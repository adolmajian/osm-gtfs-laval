import os, sys, shutil
import csv, json
import ogr, osr
import logging
from tqdm import tqdm
from collections import defaultdict, Counter
from matplotlib import pyplot as plt, figure
import copy
import numpy as np
from operator import itemgetter

import overpy

import xml.etree.ElementTree as ET

ogr.UseExceptions()
osr.UseExceptions()

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Set up logging for stdout
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter("[{asctime}][{levelname}][{name}]: {message}", style='{')
handler.setFormatter(formatter)
logging.getLogger().addHandler(handler)

api = overpy.Overpass()

relation_to_area_factor = 3600000000

region_ids = {
    "Laval": 3532125 + relation_to_area_factor,
    "Montreal": 1571328 + relation_to_area_factor
}

bus_stop_tmpl = """
    area({})->.searchArea;
    (
    node["highway"="bus_stop"](area.searchArea);
    way["highway"="platform"](area.searchArea);

    node["public_transport"="platform"]["bus"="yes"](area.searchArea);
    node["public_transport"="stop_position"]["bus"="yes"](area.searchArea);

    way["amenity"="shelter"](area.searchArea);
    node["amenity"="shelter"](area.searchArea);
    );
    out body;
"""

service_route_tmpl = """
    area({})->.searchArea;
    (
    relation["type"="route"]["route"="bus"](area.searchArea);
    );
    out body;
"""

master_route_tmpl = """
    area({})->.searchArea;
    (
    relation["type"="master_route"]["route_master"="bus"](area.searchArea);
    );
    out body;
"""


class GTFSProcessor():

    def __init__(self, gtfs_zipfile: str, boundaries_dir: str, output_dir: str):
        self.gtfs_zipfile = gtfs_zipfile
        self.output_dir = output_dir
        self.gtfs_dir = os.path.join(os.getcwd(), 'gtfs')

        # Create output directory
        logger.info('Creating output directory')
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)
        os.makedirs(self.output_dir)

        # Unzip GTFS archive
        logger.info('Unpacking GTFS archive')
        shutil.unpack_archive(self.gtfs_zipfile, self.gtfs_dir)

        # Load GTFS data into memory
        logger.info('Loading GTFS data into memory')
        self.gtfs_data = {}

        filenames = os.listdir(self.gtfs_dir)

        for filename in tqdm(filenames):
            table_name = filename[:-4]
            path = os.path.join(self.gtfs_dir, filename)
            self.gtfs_data[table_name] = {
                "path": path,
            }

            with open(path, encoding='utf-8') as csvfile:
                reader = csv.reader(csvfile)

                field_names = next(reader)
                self.gtfs_data[table_name]["field_names"] = field_names

                dict_reader = csv.DictReader(csvfile, fieldnames=field_names)
                data = [row for row in dict_reader]
                self.gtfs_data[table_name]["data"] = data

        # Load boundaries into memory
        logger.info('Loading boundary data into memory')
        self.boundaries = {}
        boundary_files = os.listdir(boundaries_dir)

        for boundary_file in tqdm(boundary_files):
            path = os.path.join(boundaries_dir, boundary_file)
            city = boundary_file[:-4]

            with open(path) as f:
                geom = ogr.CreateGeometryFromWkt(f.read())
                self.boundaries[city] = geom

    def get_latest_service_id(self):
        service_prefix = None
        logger.info('Determining latest service from calendar file')
        dates = set()
        for service in self.gtfs_data['calendar']['data']:
            dates.add(service['start_date'])

        max_date = max(list(dates))

        for service in self.gtfs_data['calendar']['data']:
            if service['start_date'] == max_date:
                service_prefix = service['service_id'][0:4]
                logger.info('Latest service is {}'.format(service_prefix))
                break

        self.service_prefix = service_prefix

    def convert_gtfs_stops_to_osm(self):
        osm_id = -100000
        gtfs_stops = []

        logger.info('Converting GTFS stops to OSM format')

        for stop in tqdm(self.gtfs_data['stops']['data']):

            if not stop['stop_id'].startswith(self.service_prefix):
                continue

            osm_id -= 1

            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(float(stop['stop_lon']), float(stop['stop_lat']))

            gtfs_stop = {
                'props': {
                    'id': osm_id,
                    'lon': stop['stop_lon'],
                    'lat': stop['stop_lat'],
                },
                'tags': {
                    'bus': 'yes',
                    'highway': 'bus_stop',
                    'name': stop['stop_name'].split('[')[0].strip(),
                    'public_transport': 'platform',
                    'ref': stop['stop_code'],
                    'shelter': 'yes' if stop['stop_abribus'] == '1' else 'no',
                },
                "gtfs_props": {
                    'stop_id': stop['stop_id'],
                    'location_type': stop['location_type'],
                    'stop_display': stop['stop_display'],
                },
                "geom": point,
            }

            gtfs_stops.append(gtfs_stop)

        logger.info('Writing converted stops to GeoJSON for visualization')
        out_path = os.path.join(output_dir, 'gtfs_stops.geojson')
        GTFSProcessor.write_data_to_geojson(gtfs_stops, out_path, "geom", ["props", "tags", "gtfs_props"])

        self.gtfs_stops = gtfs_stops

    def get_existing_osm_data(self):
        self.existing_data = {}

        existing_stops = []
        existing_routes = []
        existing_route_masters = []

        # Get existing stops
        logger.info('Getting existing stops...')
        stops_result_laval = api.query(bus_stop_tmpl.format(region_ids['Laval']))
        stops_result_montreal = api.query(bus_stop_tmpl.format(region_ids['Montreal']))

        for node in tqdm(stops_result_laval.nodes + stops_result_montreal.nodes):
            geom = ogr.Geometry(ogr.wkbPoint)
            geom.AddPoint(float(node.lon), float(node.lat))

            existing_stop = {
                "props": {
                    "id": node.id,
                    "lat": node.lat,
                    "lon": node.lon
                },
                "tags": node.tags,
                "geom": geom
            }
            existing_stops.append(existing_stop)

        # Get existing routes
        logger.info('Getting existing route relations...')
        routes_result_laval = api.query(service_route_tmpl.format(region_ids['Laval']))
        route_masters_result_laval = api.query(master_route_tmpl.format(region_ids['Laval']))

        values = ['STL', 'Société de transport de Laval']

        for relation in routes_result_laval.relations:
            tags = relation.tags
            true_count = 0

            for v in tags.values():
                if v in values:
                    true_count += 1

            if true_count > 0:

                relation_members = []

                if relation.members:
                    for member in relation.members:
                        if member._type_value == 'node':
                            relation_members.append(member)

                existing_route = {
                    "props": {
                        "id": relation.id
                    },
                    "tags": tags,
                    "members": {
                        "nodes": relation_members
                    }
                }

                existing_routes.append(existing_route)

        for relation in route_masters_result_laval.relations:
            tags = relation.tags
            true_count = 0

            for v in tags.values():
                if v in values:
                    true_count += 1

            if true_count > 0:

                relation_members = []

                if relation.members:
                    for member in relation.members:
                        if member._type_value == 'relation':
                            relation_members.append(member)

                existing_master_route = {
                    "props": {
                        "id": relation.id
                    },
                    "tags": tags,
                    "members": {
                        "relations": relation_members
                    }
                }

                existing_route_masters.append(existing_master_route)

        logger.info('Found {} existing stops in Laval and Montreal'.format(len(existing_stops)))
        logger.info('Found {} existing routes in Laval'.format(len(existing_routes)))
        logger.info('Found {} existing master routes in Laval'.format(len(existing_route_masters)))

        # Write existing stops to geojson
        logger.info('Writing existing stops to GeoJSON for visualization')
        out_path = os.path.join(output_dir, 'existing_stops.geojson')
        GTFSProcessor.write_data_to_geojson(existing_stops, out_path, "geom", ["props", "tags"])

        self.existing_data.update({
            "stops": existing_stops,
            "routes": existing_routes,
            "route_masters": existing_route_masters
        })

        logger.info('')

    def conflate_stops(self):
        # all_stops = self.gtfs_stops + self.existing_data['stops']

        final_stops = []

        logger.info('Calculating coverages for merging')
        buffer_gtfs = 5  # meter
        buffer_osm = 10  # meter
        coverage_points = ogr.Geometry(ogr.wkbMultiPoint)

        #TODO: Write extent to self for JOSM extent later

        for gtfs_stop in self.gtfs_stops:
            coverage_points.AddGeometry(gtfs_stop['geom'])

        coverage_points_utm = GTFSProcessor.reproject_geometry(coverage_points.Clone(), 4326, 32618)

        coverage_utm_gtfs = coverage_points_utm.Buffer(buffer_gtfs)
        coverage_utm_gtfs_dissolved = coverage_utm_gtfs.UnionCascaded()
        coverage_gtfs = GTFSProcessor.reproject_geometry(coverage_utm_gtfs_dissolved, 32618, 4326)

        coverage_utm_osm = coverage_points_utm.Buffer(buffer_osm)
        coverage_utm_osm_dissolved = coverage_utm_osm.UnionCascaded()
        coverage_osm = GTFSProcessor.reproject_geometry(coverage_utm_osm_dissolved, 32618, 4326)

        logger.info('Writing proximity coverage to GeoJSON')
        coverage_path_gtfs = os.path.join(self.output_dir, '{}m_coverage_gtfs.geojson'.format(buffer_gtfs))
        coverage_path_osm = os.path.join(self.output_dir, '{}m_coverage_osm.geojson'.format(buffer_osm))
        GTFSProcessor.write_geometry_to_geojson(coverage_gtfs, coverage_path_gtfs)
        GTFSProcessor.write_geometry_to_geojson(coverage_osm, coverage_path_osm)

        logger.info('Merging stops by proximity')
        for gtfs_buffer in tqdm(coverage_gtfs):

            potential_stop = None

            intersections = []
            for gtfs_stop in self.gtfs_stops:
                if gtfs_stop['geom'].Intersects(gtfs_buffer):
                    intersections.append(gtfs_stop)

            if len(intersections) > 0:
                codes = [intersection['tags']['ref'] for intersection in intersections]
                ids = [intersection['gtfs_props']['stop_id'] for intersection in intersections]

                potential_stop = intersections[0].copy()
                potential_stop['tags']['ref'] = codes
                potential_stop['gtfs_props']['stop_id'] = ids

            for osm_buffer in coverage_osm:

                if osm_buffer.Intersects(gtfs_buffer):

                    osm_intersections = []
                    for existing_stop in self.existing_data['stops']:
                        if existing_stop['geom'].Intersects(osm_buffer):
                            osm_intersections.append(existing_stop)

                    osm_int_count = len(osm_intersections)
                    if osm_int_count > 0:
                        potential_stop["props"]["id"] = osm_intersections[0]['props']['id']
                        potential_stop['props'].update({
                            "action": "modify"
                        })

                    if osm_int_count > 1:
                        logger.warning('More than one OSM stop in proximity, taking first occurence.')

            if potential_stop is None:
                logger.warning('No stop found during merging... This should not happen')
            else:
                final_stops.append(potential_stop)

        self.final_stops = final_stops

        logger.info('GTFS stops count before merging: {}'.format(len(self.gtfs_stops)))
        logger.info('GTFS stops count after merging: {}'.format(len(self.final_stops)))

        logger.info('Writing final stops to file')
        final_stops_path = os.path.join(self.output_dir, 'final_stops.geojson')
        GTFSProcessor.write_data_to_geojson(final_stops, final_stops_path, "geom", ["props", "tags", "gtfs_props"])


    def create_route_masters(self):
        osm_id_route_master = -1000
        osm_id_route = -10000

        '''
        route_master_template
        
        {
            props: {
                id: new or existing id
            }
            tags: {
                name: Henri-Bourassa - Metro Montmorency
                ref: route_short_name
                network: STL
                operator: STL
                type: route_master
                route_master: bus
                public_transport:version: 2
            }
            members: []
        }
        
        
        route_template
        {
            props: {
                id: new or existing id
            }
            tags: {
                name: Direction Montmorency
                ref: 2O
                type: route
                route: bus
                network: STL
                operator: STL
                from: name of first stop
                to: name of last stop
                round_trip: yes or no
                public_transport:version: 2
            }
            members: []
        }
        
        '''


        route_masters = defaultdict(list)
        route_master_relations = []
        route_relationss = []

        filtered_gtfs_routes = [route for route in self.gtfs_data['routes']['data'] if route['route_id'].startswith(self.service_prefix)]

        for filtered_route in filtered_gtfs_routes:
            route_masters[filtered_route['route_short_name']].append(filtered_route)

        for route_master in route_masters:
            pass



        print('')



    @staticmethod
    def write_data_to_geojson(data, out_path, geom_field, field_keys: list = None, epsg_id=None):
        # Create path
        if os.path.exists(out_path):
            os.remove(out_path)

        # Get GeoJSON driver
        driver = ogr.GetDriverByName('GeoJSON')

        ds = driver.CreateDataSource(out_path)

        spatial_ref = osr.SpatialReference()
        if epsg_id:
            spatial_ref.ImportFromEPSG(epsg_id)
        else:
            spatial_ref.ImportFromEPSG(4326)

        # Get geom type
        geom_type = data[0][geom_field].GetGeometryType()

        layer = ds.CreateLayer(out_path, geom_type=geom_type, srs=spatial_ref)

        # Create the fields
        field_names = []
        if field_keys:
            for field_key in field_keys:
                field_names.extend(data[0][field_key].keys())
            for field_name in field_names:
                layer.CreateField(ogr.FieldDefn(field_name, ogr.OFTString))  # All strings

        layer_defn = layer.GetLayerDefn()

        for item in data:
            feature = ogr.Feature(layer_defn)

            if field_keys:
                for field_key in field_keys:
                    for key in item[field_key].keys():
                        for field_name in field_names:
                            if key == field_name:
                                feature.SetField(field_name, str(item[field_key][key]))

            feature.SetGeometry(item[geom_field])

            layer.CreateFeature(feature)

    @staticmethod
    def write_geometry_to_geojson(geom, out_path):
        if os.path.exists(out_path):
            os.remove(out_path)

        driver = ogr.GetDriverByName('GeoJSON')
        ds = driver.CreateDataSource(out_path)

        geom_type = geom.GetGeometryType()

        layer = ds.CreateLayer(out_path, geom_type=geom_type)
        layer_defn = layer.GetLayerDefn()

        feature = ogr.Feature(layer_defn)
        feature.SetGeometry(geom)
        layer.CreateFeature(feature)

    @staticmethod
    def reproject_geometry(geom, in_epsg, out_epsg, return_wkt=False):
        import ogr, osr

        source = osr.SpatialReference()
        source.ImportFromEPSG(in_epsg)

        target = osr.SpatialReference()
        target.ImportFromEPSG(out_epsg)

        transform = osr.CoordinateTransformation(source, target)

        geom.Transform(transform)

        if return_wkt:
            return geom.ExportToWkt()
        else:
            return geom


gtfs_zipfile = '/home/arthur/osm-gtfs-laval/gtfs.zip'
boundaries_dir = '/home/arthur/osm-gtfs-laval/boundaries'
output_dir = '/home/arthur/osm-gtfs-laval/output'

gtfs_processor = GTFSProcessor(gtfs_zipfile, boundaries_dir, output_dir)
gtfs_processor.get_latest_service_id()
gtfs_processor.convert_gtfs_stops_to_osm()
gtfs_processor.get_existing_osm_data()
gtfs_processor.conflate_stops()
gtfs_processor.create_route_masters()

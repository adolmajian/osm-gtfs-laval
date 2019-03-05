import os, sys, shutil
import csv, json
import ogr, osr
import logging
from tqdm import tqdm
from collections import defaultdict
from operator import itemgetter
import time
import argparse

import overpy

import xml.etree.ElementTree as ET
import xml.dom.minidom as mn

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
    node["public_transport"="platform"]["bus"="yes"](area.searchArea);

    );
    out body;
"""

all_stop_tmpl = """
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
    relation["type"="route_master"]["route_master"="bus"](area.searchArea);
    );
    out body;
"""

all_tmpl = """
area({})->.searchArea;
    (
    node["highway"="bus_stop"](area.searchArea);
    relation["type"="route"]["route"="bus"]["network"="STL"](area.searchArea);
    relation["type"="route"]["route"="bus"]["network"="Société de transport de Laval"](area.searchArea);
    relation["type"="route_master"]["route_master"="bus"](area.searchArea);
    );
    out body;
"""


class GTFSProcessor():

	def __init__(self, gtfs_zipfile: str, boundaries_dir: str, output_dir: str):
		self.gtfs_zipfile = gtfs_zipfile
		self.output_dir = output_dir
		self.gtfs_dir = os.path.join(os.path.dirname(self.gtfs_zipfile), 'gtfs')

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

	def filter_gtfs_data(self):
		logger.info('Filtering the GTFS data...')
		prefix = self.service_prefix

		filtered_stops = [x for x in self.gtfs_data['stops']['data'] if x['stop_id'].startswith(prefix)]
		filtered_routes = [x for x in self.gtfs_data['routes']['data'] if x['route_id'].startswith(prefix)]
		filtered_trips = [x for x in self.gtfs_data['trips']['data'] if x['trip_id'].startswith(prefix)]
		filtered_stop_times = [x for x in self.gtfs_data['stop_times']['data'] if x['trip_id'].startswith(prefix)]

		self.gtfs_data['stops']['data'] = filtered_stops
		self.gtfs_data['routes']['data'] = filtered_routes
		self.gtfs_data['trips']['data'] = filtered_trips
		self.gtfs_data['stop_times']['data'] = filtered_stop_times

	def convert_gtfs_stops_to_osm(self):
		osm_id = -100000
		gtfs_stops = []

		logger.info('Converting GTFS stops to OSM format')

		for stop in tqdm(self.gtfs_data['stops']['data']):
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
		time.sleep(30)
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

				relation_member_nodes = []
				relation_member_ways = []

				if relation.members:
					for member in relation.members:
						if member._type_value == 'node':
							relation_member_nodes.append({
								"props": {
									"type": 'node',
									"ref": member.ref,
									"role": member.role
								}
							})
						if member._type_value == 'way':
							relation_member_ways.append({
								"props": {
									"type": 'way',
									"ref": member.ref,
									"role": ""
								}
							})

				existing_route = {
					"props": {
						"id": relation.id
					},
					"tags": tags,
					"members": {
						"nodes": relation_member_nodes,
						"ways": relation_member_ways
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

				relation_member_nodes = []

				if relation.members:
					for member in relation.members:
						if member._type_value == 'relation':
							relation_member_nodes.append(member)

				existing_master_route = {
					"props": {
						"id": relation.id
					},
					"tags": tags,
					"members": {
						"relations": relation_member_nodes
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

		logger.info('... {} existings OSM stops found'.format(len(existing_stops)))
		logger.info('... {} existings OSM routes found'.format(len(existing_routes)))
		logger.info('... {} existings OSM route masters found'.format(len(existing_route_masters)))

		self.existing_data.update({
			"stops": existing_stops,
			"routes": existing_routes,
			"route_masters": existing_route_masters
		})

	def write_route_ids_csv(self):
		existing_routes = self.existing_data['routes']

		def complete(name):
			vals = ['N', 'S', 'E', 'O']

			b = False

			for val in vals:
				if val in name:
					b = True
					break
			return b

		rows = []
		for existing_route in existing_routes:
			osm_id = existing_route['props']['id']
			ref = existing_route['tags']['ref']
			name = int(ref[:-1]) if complete(ref) else -99999
			rows.append([osm_id, name, ref])

		rows.sort(key=lambda x: x[1])

		logger.info('Writing existing route ids to CSV')
		existing_dir = os.path.join(self.output_dir, 'routes')
		os.makedirs(existing_dir)
		out_path = os.path.join(existing_dir, 'existing_routes.csv')

		with open(out_path, 'w') as f:
			writer = csv.writer(f)
			writer.writerow(['osm_id', 'route', 'ref'])

			for row in rows:
				writer.writerow(row)

	def conflate_stops(self):
		final_stops = []

		logger.info('Calculating coverages for merging')
		buffer_gtfs = 5  # meter
		buffer_osm = 5  # meter
		coverage_points = ogr.Geometry(ogr.wkbMultiPoint)

		for gtfs_stop in self.gtfs_stops:
			coverage_points.AddGeometry(gtfs_stop['geom'])

		self.extent = coverage_points.GetEnvelope()

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

		logger.info('Filtering coverages for Laval')
		coverage_osm_laval = []
		for osm_buffer in coverage_osm:
			if osm_buffer.Intersects(self.boundaries['Laval']):
				coverage_osm_laval.append(osm_buffer)

		logger.info('Merging stops by proximity')
		for gtfs_buffer in tqdm(coverage_gtfs):

			potential_stop = None

			intersections = []
			for gtfs_stop in self.gtfs_stops:
				if gtfs_stop['geom'].Intersects(gtfs_buffer):
					intersections.append(gtfs_stop)

			if len(intersections) > 0:
				codes = list(set([intersection['tags']['ref'] for intersection in intersections]))
				codes_joined = ";".join(codes)
				ids = list(set([intersection['gtfs_props']['stop_id'] for intersection in intersections]))

				potential_stop = intersections[0].copy()
				potential_stop['tags']['ref'] = codes_joined
				potential_stop['gtfs_props']['stop_id'] = ids

			for osm_buffer in coverage_osm_laval:

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
						logger.warning('More than one OSM stop in proximity, taking first occurrence.')
						print(osm_intersections)

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

	def create_relations(self):
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
				ref: 2E
				type: route
				route: bus
				network: STL
				operator: STL
				from: name of first stop
				to: name of last stop
				roundtrip: yes or no
				public_transport:version: 2
			}
			members: []
		}
		
		'''

		route_masters = defaultdict(list)
		route_master_relations = []
		route_relations = []

		# Create directory for shapes
		shapes_dir = os.path.join(self.output_dir, 'shapes')
		if os.path.exists(shapes_dir):
			shutil.rmtree(shapes_dir)
		os.makedirs(shapes_dir)

		gtfs_routes = self.gtfs_data['routes']['data']

		for gtfs_route in gtfs_routes:
			route_masters[gtfs_route['route_short_name']].append(gtfs_route)

		for key, route_master in tqdm(route_masters.items()):
			logger.info('Creating route relations for {}'.format(key))
			round_trip = "yes" if len(route_master) == 1 else "no"

			member_routes = []

			for route in route_master:
				osm_id_route -= 1

				route_ref = route['route_url'].split('route_id=')[1].strip()
				route_name = route['route_long_name']

				member_nodes = []
				# member_ways = []

				# logger.info('... Finding longest trip and stops for {}'.format(route_ref))
				first_stop_name, last_stop_name, stop_ids, trip_id = self.get_route_stops(route['route_id'])

				# Build the shape to write to geojson
				trip_shape = self.make_trip_shape(trip_id)
				trip_shape.update({
					"fields": {
						"ref": route_ref
					}
				})
				trip_filename = os.path.join(shapes_dir, '{}_route.geojson'.format(route_ref))
				GTFSProcessor.write_data_to_geojson([trip_shape], trip_filename, 'geom', ['fields'])

				for i, stop_id in enumerate(stop_ids):
					member_node = {
						"props": {
							"type": 'node',
							"ref": stop_id,
							"role": 'platform'
						}
					}
					# if i == 0:
					# 	member["props"]["role"] = 'platform_entry_only'
					# if i == len(stop_ids) - 1:
					# 	member['props']['role'] = 'platform_exit_only'
					member_nodes.append(member_node)

				route_relation = {
					"props": {
						"id": osm_id_route
					},
					"tags": {
						"name": route_name,
						"ref": route_ref,
						"type": "route",
						"route": "bus",
						"network": "STL",
						"operator": "STL",
						"from": first_stop_name,
						"to": last_stop_name,
						"roundtrip": round_trip,
						"public_transport:version": 2
					},
					"members": {
						"nodes": member_nodes
					}
				}
				route_relations.append(route_relation)

				member_route = {
					"props": {
						"type": "relation",
						"ref": osm_id_route,
						"role": ""
					}
				}
				member_routes.append(member_route)

			# Create route master relation
			# logger.info('Creating route master relation for {}'.format(key))

			# Only create route_master_relation if there are more than 1 directions
			if round_trip == "yes":
				logger.info('... {} is probably a loop route. Skipping.'.format(key))
				continue

			osm_id_route_master -= 1

			directions = [x['route_long_name'].split('Direction')[1].strip() for x in route_master]
			name = " - ".join(directions)

			route_master_relation = {
				"props": {
					"id": osm_id_route_master
				},
				"tags": {
					"name": name,
					"ref": key,
					"network": "STL",
					"operator": "STL",
					"type": "route_master",
					"route_master": "bus",
					"public_transport:version": 2
				},
				"members": member_routes
			}
			route_master_relations.append(route_master_relation)

			self.route_relations = route_relations
			self.route_master_relations = route_master_relations

	def conflate_relations(self):
		logger.info('Resolving route relation conflicts')
		for existing_route in self.existing_data['routes']:
			existing_ref = existing_route['tags']['ref']
			existing_name = existing_route['tags']['name']

			for route_relation in self.route_relations:
				new_ref = route_relation['tags']['ref']
				new_name = route_relation['tags']['name']
				new_id = route_relation["props"]["id"]

				if new_ref == existing_ref and new_name == existing_name:
					logger.info('... Found existing route. Resolving.')
					# Set id to existing id
					# add action=modify to props
					# Replace the tags with new tags
					# Replace node relation list with new ones
					# Loop through master and replace ref values with existing id
					# Append way members if any

					existing_id = existing_route["props"]["id"]

					route_relation["props"]["id"] = existing_id
					route_relation["props"]["action"] = "modify"

					if 'ways' in existing_route['members']:
						route_relation['members']['ways'] = existing_route['members']['ways']

					for route_master_relation in self.route_master_relations:
						if route_master_relation['tags']['ref'] == new_ref[:-1]:
							for member in route_master_relation['members']:
								if member['props']['ref'] == new_id:
									member['props']['ref'] = existing_id
									break
					break

	def write_to_xml(self, types: list):
		print(types)
		logger.info('Writing JOSM XML files for {}.'.format(', '.join(types)))
		output_file = os.path.join(self.output_dir, 'gtfs_laval.xml')

		if os.path.exists(output_file):
			os.remove(output_file)

		root = ET.Element("osm")
		root.set('version', '0.6')

		min_lon, max_lon, min_lat, max_lat = self.extent

		bounds = ET.SubElement(root, 'bounds')
		bounds.set('minlat', str(min_lat))
		bounds.set('minlon', str(min_lon))
		bounds.set('maxlat', str(max_lat))
		bounds.set('maxlon', str(max_lon))

		if 'stops' in types:
			for stop in self.final_stops:
				node = ET.SubElement(root, 'node')
				node.set('version', '1')

				for k, v in stop['tags'].items():
					tag = ET.SubElement(node, 'tag')
					tag.set('k', str(k))
					tag.set('v', str(v))

				for k, v in stop['props'].items():
					node.set(k, str(v))

		if 'routes' in types:
			for route_relation in self.route_relations:
				relation = ET.SubElement(root, 'relation')
				relation.set('version', '1')

				for k, v in route_relation['tags'].items():
					tag = ET.SubElement(relation, 'tag')
					tag.set('k', str(k))
					tag.set('v', str(v))

				for k, v in route_relation['props'].items():
					relation.set(k, str(v))

				member_types = route_relation['members'].keys()

				for member_node in route_relation['members']['nodes']:
					mem = ET.SubElement(relation, 'member')

					for k, v in member_node['props'].items():
						mem.set(k, str(v))

				if 'ways' in member_types:
					for member_way in route_relation['members']['ways']:
						mem = ET.SubElement(relation, 'member')

						for k, v in member_way['props'].items():
							mem.set(k, str(v))

		if 'route_masters' in types:
			for route_master_relation in self.route_master_relations:
				relation = ET.SubElement(root, 'relation')
				relation.set('version', '1')

				for key, value in route_master_relation.items():
					if key == 'tags':
						for k, v in route_master_relation['tags'].items():
							tag = ET.SubElement(relation, 'tag')
							tag.set('k', str(k))
							tag.set('v', str(v))
					if key == 'props':
						for k, v in route_master_relation['props'].items():
							relation.set(k, str(v))
					if key == 'members':
						for member in route_master_relation['members']:
							mem = ET.SubElement(relation, 'member')

							for key, value in member.items():
								if key == 'props':
									for k, v in member['props'].items():
										mem.set(k, str(v))

		tree = ET.ElementTree(root)
		tree.write(output_file, encoding='unicode')

	def get_route_stops(self, route_id):
		trips = defaultdict(dict)

		for trip in self.gtfs_data['trips']['data']:
			if trip['route_id'] == route_id:
				trip_id = trip['trip_id']

				stops = []

				for stop_time in self.gtfs_data['stop_times']['data']:
					if stop_time['trip_id'] == trip_id:
						stop = {
							"stop_id": stop_time["stop_id"],
							"stop_sequence": int(stop_time["stop_sequence"])
						}

						stops.append(stop)

				trips[trip_id]['stops'] = stops
				trips[trip_id]['stop_count'] = len(stops)

		trips = dict(trips)
		# logger.info('... ... Found {} trips for route id {}'.format(len(trips.keys()), route_id))

		# # Find the longest trip and sort the stops
		longest_trip_id = max(trips, key=lambda v: trips[v]['stop_count'])
		# logger.info('... ... ... the longest trip is {} with {} stops'.format(longest_trip_id,
		# 																	  trips[longest_trip_id]['stop_count']))

		# Sort the stops
		stops = trips[longest_trip_id]['stops']
		# stops.sort(key=itemgetter('stop_sequence'))
		stops_sorted = sorted(stops, key=itemgetter('stop_sequence'))

		# Map to final stop ids
		first_stop_name = None
		last_stop_name = None
		stop_ids = []
		for i, stop in enumerate(stops_sorted):
			for final_stop in self.final_stops:
				if stop['stop_id'] in final_stop['gtfs_props']['stop_id']:
					stop_ids.append(final_stop['props']['id'])
					if i == 0:
						first_stop_name = final_stop['tags']['name']
					if i == len(stops) - 1:
						last_stop_name = final_stop['tags']['name']

		return first_stop_name, last_stop_name, stop_ids, longest_trip_id

	def make_trip_shape(self, trip_id):
		trips = self.gtfs_data['trips']['data']
		shapes = self.gtfs_data['shapes']['data']

		shape_id = None
		for trip in trips:
			if trip['trip_id'] == trip_id:
				shape_id = trip['shape_id']

		shape_points = []
		for shape in shapes:
			if shape['shape_id'] == shape_id:
				shape_point = {
					"lon": float(shape['shape_pt_lon']),
					"lat": float(shape['shape_pt_lat']),
					"sequence": int(shape['shape_pt_sequence'])
				}
				shape_points.append(shape_point)

		shape_points.sort(key=itemgetter('sequence'))

		line = ogr.Geometry(ogr.wkbLineString)
		for shape_point in shape_points:
			line.AddPoint(shape_point['lon'], shape_point['lat'])

		return {"geom": line}

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

cwd = os.getcwd()
gtfs_zipfile = os.path.join(cwd, 'gtfs.zip')
boundaries_dir = os.path.join(cwd, 'boundaries')
output_dir = os.path.join(cwd, 'output')

time1 = time.time()
gtfs_processor = GTFSProcessor(gtfs_zipfile, boundaries_dir, output_dir)
gtfs_processor.get_latest_service_id()
gtfs_processor.filter_gtfs_data()
gtfs_processor.convert_gtfs_stops_to_osm()
gtfs_processor.get_existing_osm_data()
gtfs_processor.write_route_ids_csv()
gtfs_processor.conflate_stops()
gtfs_processor.create_relations()
gtfs_processor.conflate_relations()
gtfs_processor.write_to_xml(['route_masters'])

time2 = time.time()
duration = time2 - time1
logger.info('Processing completed in {:.2f} seconds'.format(duration))

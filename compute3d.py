# coding: utf-8
import sys
sys.path.append('.')
import os

import json

from pydmga.geometry import OrthogonalGeometry
from pydmga.container import Container
from pydmga.diagram import Diagram
from pydmga.io import pdb
from pydmga.io import helper
import traceback

from pydmga.draw import assets
from pydmga.draw import render
from pydmga.draw.render import Asset
from pydmga.draw.render import Color

from math import pi
from math import sqrt
from math import cos
from math import sin
from math import acos
from math import asin

import numpy as np

###############################################################
# REAL SCRIPT STARTS HERE #####################################
###############################################################

# TODO: Dipol angle (lists per time + sensowny histogram)
# TODO: Save layer as a .ndx file (each frame as a single group?), to be able to recover layers
# TODO: save neighbours list in TXT readable file for later analysis and/or speed up?
# TODO: rethink computation of following layers by providing .ndx file (one entry per frame), so that
#       we can use .ndx files generated by one run for base molecules to compute next layers
#       in case of a single group in .ndx use it for all frames (handy for base molecules around which water flows)

# # IMPORT SETTINGS
# # WE WILL USE THEM IN FUNCTIONS?
if (len(sys.argv) < 2):
	print "Computes the first hydratation layer statistics."
	print ""
	print "Usage: python compute3d.py [path.to.conf] {no-backup}"
	print ""
	print "  [path.to.conf] is a python module path to a file containing configuration.";
	print "                 see sample files in jobs/* for reference."
	print "  {no-backup}    optional. If theres is no-backup in the end, "
	print "                 then the output directory will not be backuped."
	print "                 By default for safety, we always do backup."
	print "                 No backup is good for recomputing or computing to "
	print "                 common directory by many scripts."
	sys.exit(-1)

job = helper.Job(sys.argv, app="compute3d")
outf_basename = job.settings.OUTPUT_DIRPATH

# common job.settings.
OTHER_RADIUS = 1.0
max_analys_frames = 100
max_voro_frames = 10
molecules = {
	'chol': {
		'def': {
			'lo': 1,
			'up': 74,
			'plane_atoms_ids': [9, 18, 35],
			# we fit consecutive MOLECULES of this type to this atoms 
			# we want:
			# a) all atoms on the same plane
			# b) line from plane_atoms_ids[1] to plane_atoms_ids[2] on time t-1 be on the line at t )
			# c) position of plane_atoms_ids[2] on t-1 the same as in t)
		},
	},
}
try: max_analys_frames = job.settings.MAX_ANALYSIS_FRAMES
except: job.log("# WARNING: using default MAX_ANALYSIS_FRAMES = {}\n".format(max_analys_frames))
try: max_voro_frames = job.settings.MAX_VORO_FRAMES 
except: job.log("# WARNING: using default MAX_VORO_FRAMES = {}\n".format(max_voro_frames))
try: molecules = job.settings.MOLECULES; 
except: job.log("# WARNING: using default MOLECULES (CHOL) = {}\n".format(str(molecules)))

for key, moldata in molecules.iteritems():
	if "ndx" not in moldata:
		if "def" in moldata:
			moldata["ndx"] = list(range(moldata['def']["lo"], moldata['def']["up"] + 1))
		else:
			raise Exception("You must specify either 'def' or 'ndx' for '{}' in settings.".format(key))

# make a dict that translates ATOM_ID to MOLECULE key. 
id_to_mol = dict()
for key, moldata in molecules.iteritems():
	for i in range(moldata['def']['lo'], moldata['def']['up']+1):
		id_to_mol[i] = key


def id_to_molecule(i):
	if i in id_to_mol:
		return id_to_mol[i]
	return None


def num_to_molecule(i):
	if (i+1) in id_to_mol:
		return id_to_mol[i+1]
	return None


def is_hydrogen(item):
	# TODO: THIS IS TOO SPECYFIC!
	return item["atom"][0] == 'H' or item["atom"][1] == 'H'


in_file = file(job.settings.PDB_TRAJ_FILEPATH, "r")  # read from here

# generate as many analysis files as there are different kinds of molecules in the setup
for molecule_key in molecules:	
	full_graph_file = job.user_file("{molecule_key}_{app}_full_graph.tjson", {"molecule_key": molecule_key})
	molecules[molecule_key]['full_graph_file'] = full_graph_file
	water_neighbours_file = job.user_file("{molecule_key}_{app}_water_neighbours.tjson", {"molecule_key": molecule_key})
	molecules[molecule_key]['water_neighbours_file'] = water_neighbours_file
	ion_neighbours_file = job.user_file("{molecule_key}_{app}_ion_neighbours.tjson", {"molecule_key": molecule_key})
	molecules[molecule_key]['ion_neighbours_file'] = ion_neighbours_file	
	internal_molecule_neighbours_file = job.user_file("{molecule_key}_{app}_internal_molecule_neighbours.tjson", {"molecule_key": molecule_key})
	external_molecule_neighbours_file = job.user_file("{molecule_key}_{app}_external_molecule_neighbours.tjson", {"molecule_key": molecule_key})
	molecules[molecule_key]['internal_molecule_neighbours_file'] = internal_molecule_neighbours_file
	molecules[molecule_key]['external_molecule_neighbours_file'] = external_molecule_neighbours_file

	for group in ['external_molecule_neighbours', 'ion_neighbours', 'water_neighbours']:
		# write here analysis data per frame
		group_file = job.user_file("{molecule_key}_{app}_contact_{group}.dat", {"molecule_key": molecule_key, "group": group})
		group_file.write(job.info)
		group_file.write("# LEGEND: \n")
		group_file.write("# time frame              obvious\n")
		group_file.write("# avg/max/min/sum[dist]   apropriate function of a computed distance from molecule atom\n")
		group_file.write("#                         to group molecules neighbours in Voronoi diagram. \n")
		group_file.write("#                         this is half the distance between centers of the atoms (in the basic application).\n")
		group_file.write("#                         from atom to voronoi cell face.\n")
		group_file.write("# area                    sum of the contact area with group.\n")
		group_file.write("# \n")
		group_file.write("# This file was generated for job '{}', molecule: '{}', and group '{}'\n".format(job.jobname, molecule_key, group))
		group_file.write("# \n")
		group_file.write("# Header numbers for help in gnuplot! \n")
		group_file.write("# 1    2     3       4       5       6       7   \n")
		group_file.write("# time frame avgdist maxdist mindist sumdist area\n")		
		molecules[molecule_key]['contact_analys_{}'.format(group)] = group_file

	volume_file = job.user_file("{molecule_key}_{app}_volume_analys.dat", {"molecule_key": molecule_key})
	volume_file.write(job.info)
	volume_file.write("# LEGEND: \n")
	volume_file.write("# time frame              obvious\n")
	volume_file.write("# volume                  sum of the voronoi volumes associated to atoms of the molecule\n")
	volume_file.write("# \n")
	volume_file.write("# This file was generated for job '{}' and molecule: '{}'\n".format(job.jobname, molecule_key))
	volume_file.write("# \n")
	volume_file.write("# Header numbers for help in gnuplot! \n")
	volume_file.write("# 1    2     3     \n")
	volume_file.write("# time frame volume\n")	
	molecules[molecule_key]['volume_analys'] = volume_file

	# TODO: in more generla version of this script make groups defined by the user, ie.
	# TODO: user provides ndx files (each group is for a given time) and the graph is between those groups (defined by user)

water_molecules_file = job.user_file("{app}_water_molecules.json")
positions_file = job.user_file("{app}_positions.tjson")

def graph_data(g, pos): 
	'''
	assumes g graph of the shape:
	{1: {5: (9, 8)}, 2: {3: (1, 2), 4: (3, 4)}}
	e.g. vertex 1 is connected to 5, with data (9, 8)
	it computes sum of all data at pos from all the neighbours
	e.g. sum_graph_data(g, 0) = [9, 1, 3]
	'''
	return sum(map(lambda x: [nbdata[pos] for nb, nbdata in x.iteritems()], [vdata for v, vdata in g.iteritems()]), [])

try:
	frame = 0
	time = 0

	water_molecules = {}
	current_water = None
	water_ready = False;

	hao_id_to_num = {}
	curr_num = 0
	hao_ready = False	

	while True: # iterate over all frames	

		# THIS IS AN OLD SCRIPT SO IT NEED TO PARSE PDB BY ITS OWN...
		# PARSING STARTS HERE...
		time = float(frame)

		while True:
			line = in_file.readline()        
			if not line:
				raise StopIteration

			dat = pdb.tokenize(line)
			if dat is None: 
				continue;

			item = pdb.PDBRecord(dat)
			if item.type() == "TITLE":
				if item["time"]:
					time = item["time"]
					print "Found frame time t =", time

			if item.type() == "CRYST1":	# create new geometry
				(a, b, c) = (item["a"], item["b"], item["c"])
				geometry = OrthogonalGeometry(item["a"], item["b"], item["c"], True, True, True)				
				container = Container(geometry)				
				print "Geometry and Container created... "				
			elif item.type() == "ATOM":			
				if (not water_ready) and job.settings.is_water(item):
					# we assume oxygens are first, then two hydrogens in the file
					if is_hydrogen(item): 
						current_water.append(item.as_particle())
					else:
						# is oxygen then
						current_water = [item.as_particle(), ]
						water_molecules[item["id"]] = current_water
				# we dump hydrogen atoms, they are not taken into account in current script! (setting hao is not neccessary now)
				if not is_hydrogen(item):
					(atom_id, x, y, z, r) = item.as_particle()
					if not hao_ready:
						hao_id_to_num[atom_id] = curr_num
						curr_num += 1
					container.add(atom_id, x, y, z, OTHER_RADIUS) # for now fix r (carbon and others)
			elif item.type() == "TER":				
				print "TER from line", line.strip(), ", frame", frame, " is ready."

				if not water_ready:
					# I save only water connections of oxy to hydrogen, as positions will be elsewhere
					# This is to save space on the disk.
					only_connections = dict(
						map(
							# this is ugly, but iteritems returns pairs of key, list
							# list has 3 entries: O, H1, H2, each entry has id in the [0] entry.
							lambda item: (item[0], (item[1][1][0], item[1][2][0])), 
							water_molecules.iteritems()
						)
					)
					json.dump(only_connections, water_molecules_file)
					water_molecules_file.close()				
					water_ready = True

				hao_ready = True
				break
			elif item.type() == "ENDMDL":
				pass
			else:
				pass
		# OK, NOW WE HAVE PARSED NEXT PDB FRAME, we can do computations.		

		# tu sie buduje Voronoi
		diagram = Diagram(container)

		# setup items to be update later
		for key, moldata in molecules.iteritems():
			# since we are inside a loop over molecules, I can initialize important variables
			# that are local to this frame diagram and individual for each molecule
			moldata['distances'] = []
			moldata['volumes'] = []
			moldata['water_neighbours'] = {}
			moldata['ion_neighbours'] = {}
			moldata['external_molecule_neighbours'] = {}
			moldata['internal_molecule_neighbours'] = {}
			moldata['full_graph'] = {}

		# in positions, I will write positions of important molecules (e.g. neighbours)
		# later, it would be better to read trajectory alongside the network file
		# and get positions from there... (TODO: jak juz bede mia?? obs??ug?? ndx / xtc / trr)
		# positions are common for all the molecules
		positions = {}

		for molkey, moldata in molecules.iteritems():
			for atom_id in moldata["ndx"]:
				if atom_id not in hao_id_to_num:
					# atom_id is not Heavy Atom (so is hydrogen) - we skip those.
					continue
				atom_num = hao_id_to_num[atom_id]
				cell = diagram.get_cell(atom_num)				
				this = container.get(atom_num)
				assert this[0] == atom_id, "{} != {}".format(this[0], atom_id)

				positions[atom_id] = tuple(this[1:4])
				p_this = np.array(this[1:4])

				# update volue of a given molecule
				moldata['volumes'].append( (atom_id, cell.volume()) )
				full_neighbours = moldata['full_graph'].setdefault(atom_id, {})								
				
				# iterate over neighbours by their cells				
				for side in cell.sides:
					other = container.get(side.neighbour)
					neighbour_id = other[0]	
					positions[neighbour_id] = tuple(other[1:4])				
					p_other = np.array(other[1:4])
					# calculate distance										
					n = float(np.linalg.norm(p_this - p_other))
					# Liczy odleglosc od ??ciany Voronoia. Jesli uzywamy tylko atomow o 
					# tym samym radiusie to dwie liczby this i other si?? zcanceluj??.
					dist = 0.5 * (n + this[4]**2 - other[4]**2)  
					# compute contact area
					area = side.area()

					# I try for the graph to be minimalistic as possible!
					# I add this side neighbour as a neighbour of the molecule graph
					full_neighbours[neighbour_id] = (dist, area)

					# now, we differentiate the types of neighbours
					if job.settings.is_water(neighbour_id): # tylko woda jako sasiad					
						neighbour = moldata['water_neighbours'].setdefault(neighbour_id, {})				
						# neighbour is a collection of neighbours in 'this' molecule of the 'other' water
						neighbour[atom_id] = (dist, area)

						# add water neighbours for the eventual water-molecule contact interaction calculation
						O, H1, H2 = water_molecules[neighbour_id]
						positions.update({
							O[0]: tuple(O[1:4]),
							H1[0]: tuple(H1[1:4]),
							H2[0]: tuple(H2[1:4]),
						})
					elif job.settings.is_ion(neighbour_id): # some ion as a neighbour
						neighbour = moldata['ion_neighbours'].setdefault(neighbour_id, {})				
						neighbour[atom_id] = (dist, area)
					else: # atom of this or the other molecule in the set					
						other_molkey = id_to_molecule(other[0])
						if other_molkey:
							if other_molkey == molkey:
								neighbour = moldata['internal_molecule_neighbours'].setdefault(neighbour_id, {})				
								neighbour[atom_id] = (dist, area)
							else:
								neighbour = moldata['external_molecule_neighbours'].setdefault(neighbour_id, {})				
								neighbour[atom_id] = (dist, area, other_molkey)
						else:
							# miscelanous molecule contact...
							pass

					# now we go to the next cell

				# next atom in this molecule

			# now we have parsed all atoms in this molecule, we can save partial results
			moldata["full_graph_file"].write("{} {}\n".format(time, json.dumps(moldata["full_graph"])))	
			moldata["ion_neighbours_file"].write("{} {}\n".format(time, json.dumps(moldata["ion_neighbours"])))	
			moldata["water_neighbours_file"].write("{} {}\n".format(time, json.dumps(moldata["water_neighbours"])))	
			moldata["internal_molecule_neighbours_file"].write("{} {}\n".format(time, json.dumps(moldata["internal_molecule_neighbours"])))	
			moldata["external_molecule_neighbours_file"].write("{} {}\n".format(time, json.dumps(moldata["external_molecule_neighbours"])))	
			# next molecule	

		# save data before we go to next frame
		positions_file.write("{} {}\n".format(time, json.dumps(positions)))			

		# output some basic statistics
		for key, moldata in molecules.iteritems():
			volume = sum(map(lambda x: x[1], moldata["volumes"]), 0.0)
			moldata['volume_analys'].write("{} {} {}\n".format(time, frame, volume))						
			for group in ['ion_neighbours', 'water_neighbours', 'external_molecule_neighbours']:
				# time, frame, avgdist, maxdist, mindist, sumdist, volume, area
				group_distances = graph_data(moldata[group], 0)
				group_areas = graph_data(moldata[group], 1)
				moldata['contact_analys_{}'.format(group)].write(
					"{} {} {} {} {} {} {}\n".format(
						time, frame, 
						sum(group_distances or [0.0]) / (len(group_distances) or 1.0), # avg
						max(group_distances or [0.0]), # max
						min(group_distances or [0.0]), # min
						sum(group_distances or [0.0]), # sum
						sum(group_areas or [0.0]),
					)
				)			

		frame += 1

except Exception as e:
	job.finish("By exception")
	print e
	traceback.print_exc()
else:
	job.finish("Finished normally")

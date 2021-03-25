# coding: utf-8
import sys
sys.path.append('.')
import os

# TODO: TU SKONCZYLEM Z TYM PLIKIEM :
# - dodac wypisywanie pliku .pdb z molekulami i ich wspolnymi scianami? (daloby sie, mysle prosto...)
# - pokomentowac lepiej fragementy...
# - to jest wielki nieczytelny skrypt, przemyslec i przepisac z dobrymi praktykami programistycznymi
#   jezeli bedziemy uzywac tego w przyszlosci

from pydmga.geometry import OrthogonalGeometry
from pydmga.container import Container
from pydmga.diagram import Diagram
from pydmga.io import pdb
import traceback

from math import pi
from math import sqrt
from math import cos
from math import sin
from math import acos
from math import asin
from collections import deque

import numpy as np	
import json

# IMPORT SETTINGS
import importlib
import datetime
if len(sys.argv) < 2:
	print "Usage: python layers.py [path.to.conf]"
	print ""
	print "  [path.to.conf] is a python module path to a file containing configuration.";
	print "                 see sample files in jobs/* for reference."
	print "  {backup_mode}  default: backup; possible values: backup|clear|override"
	print "                 backup   - make backup of previous OUTPUT_DIR"   
	print "                 override - do nothing with OUTPUT_DIR, overrride files"  	
	print ""
	sys.exit(-1)

jobname = sys.argv[1]
settings = importlib.import_module(sys.argv[1])
jobstart = datetime.datetime.now()

backup_mode = "backup"
if len(sys.argv) >= 3:
	backup_mode = sys.argv[2]
	if backup_mode not in ["backup", "override", "no-backup"]:
		backup_mode = "backup"
# END IMPORT SETTINGS.
if backup_mode == "no-backup":
	backup_mode = "override"

###############################################################
# REAL SCRIPT STARTS HERE #####################################
###############################################################
if not(settings.OUTPUT_DIRPATH):
	raise Exception("No OUTPUT_DIRPATH set in settings. BREAK.")

if os.path.exists(settings.OUTPUT_DIRPATH):
	if backup_mode == "backup":	
		outf_backup = settings.OUTPUT_DIRPATH
		counter = 0;
		while os.path.exists(outf_backup):	
			outf_backup = "{0}_{1}.bak".format(settings.OUTPUT_DIRPATH, counter)
			counter += 1
		os.rename(settings.OUTPUT_DIRPATH, outf_backup)
		os.mkdir(settings.OUTPUT_DIRPATH)
else:
	os.mkdir(settings.OUTPUT_DIRPATH)

# I use outf_basename bo chce miec wsteczna kompatybilnosc
outf_basename = settings.OUTPUT_DIRPATH

log_filepath = os.path.join(settings.OUTPUT_DIRPATH, "log_hydbond_layers.txt")
log_file = file(log_filepath, "w")


def mkdir_p(dirpath):
	try:
		os.mkdir(dirpath)
	except Exception as e:
		pass


LAYERS_DIRPATH = os.path.join(settings.OUTPUT_DIRPATH, "layers_" + jobname)
mkdir_p(LAYERS_DIRPATH)

# common settings:
HYDROGEN_RADIUS = 1.09
OTHER_RADIUS = 1.7
solv_hao = True
mols_hao = True	
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
try: solv_hao = settings.SOLVENT_HAO
except: log_file.write("# WARNING: using default SOLVENT_HAO = {}\n".format(solv_hao))
try: mols_hao = settings.MOLECULES_HAO
except: log_file.write("# WARNING: using default MOLECULES_HAO = {}\n".format(mols_hao))
try: max_analys_frames = settings.MAX_ANALYSIS_FRAMES
except: log_file.write("# WARNING: using default MAX_ANALYSIS_FRAMES = {}\n".format(max_analys_frames))
try: max_voro_frames = settings.MAX_VORO_FRAMES 
except: log_file.write("# WARNING: using default MAX_VORO_FRAMES = {}\n".format(max_voro_frames))
try: molecules = settings.MOLECULES; 
except: log_file.write("# WARNING: using default MOLECULES (CHOL) = {}\n".format(str(molecules)))
COLLECT_LAYERS = [1, 2, 3]
try: COLLECT_LAYERS = settings.COLLECT_LAYERS; 
except: log_file.write("# WARNING: using default COLLECT_LAYERS = {}\n".format(str(COLLECT_LAYERS)))
CROSS_LAYERS = [2, 3, 4]
try: CROSS_LAYERS = settings.CROSS_LAYERS; 
except: log_file.write("# WARNING: using default CROSS_LAYERS = {}\n".format(str(CROSS_LAYERS)))
COMMON_CROSS_LAYERS = [len(molecules) + i for i in range(3)]
try: COMMON_CROSS_LAYERS = settings.COMMON_CROSS_LAYERS; 
except: log_file.write("# WARNING: using default COMMON_CROSS_LAYERS = {}\n".format(str(COMMON_CROSS_LAYERS)))

#MINUS_INFTY = -((float(max(CROSS_LAYERS))+2.0)**2)
#MINUS_INFTY = - (max(CROSS_LAYERS)+2)**2
MINUS_INFTY = -1000


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
	return item["atom"][0] == 'H' or item["atom"][1] == 'H'


def to_coords(atom_item):
	# atom_item[0] is ID
	return [atom_item[1], atom_item[2], atom_item[3]]


def dist2(x, y):
	return sum(map(lambda z: (z[0]-z[1])**2, zip(x, y)))


def dist(x, y):
	return sqrt(x, y)


def pbc_positions(x, box):
    return [
        [x[0] + d[0], x[1] + d[1], x[2] + d[2]]
        for d in [
            (a, b, c)
            for a in (-box[0], 0, +box[0])
            for b in (-box[1], 0, +box[1])
            for c in (-box[2], 0, +box[2])
        ]
    ]


def closest_distance2(x, y, box):
	pbcs = pbc_positions(y, box)
	v = min([dist2(x, p) for p in pbcs])
	return v


def closest_distance(x, y, box):
	return sqrt(closest_distance2(x,y,box))


def closest_distance2_copy(x, y, box):
	pbcs = pbc_positions(y, box)
	v = min([(dist2(x, p), p) for p in pbcs])
	return v


def closest_distance_copy(x, y, box):
	d, pos = closest_distance2_copy(x,y,box)
	return sqrt(d), pos


def run_BFS(n, diagram, initials, depth, volumes, nonwater_neighbours, water_neighbours):
	"""
	Run a breadth first search starting at given set of atoms
	The connection in graph is defined per voronoi face, i.e. 
	two atoms are connected in the graph if they share a common 
	voronoi face.
	"""
	q = deque(initials)
	visited = [MINUS_INFTY for i in range(n)]
	for v in initials:
		visited[v] = 0
	while len(q) > 0:
		v = q.popleft()	
		cell = diagram.get_cell(v)		
		if volumes[v] < 0.0:
			volumes[v] = cell.volume()
		if visited[v] < depth:
			for s in cell.sides:					
				u = s.neighbour
				other = container.get(u)
				neighbour_id = other[0]
				if settings.is_water(neighbour_id): 
					if visited[u] < 0:
						visited[u] = visited[v] + 1
						q.append(u)
					if visited[v] == 0:
						# for the initials, we add to water surface, 
						# we always add, as we need to revisit 
						# all sides that touch water and single water tan touch many atoms in the initials
						# (but one water and wone atom = only one possible common side - see definition of Voronoi)
						water_neighbours.append( (v, u, neighbour_id, id_to_molecule(neighbour_id), s.area()) )
				else:
					# not water neighbour
					# we measure only the distance from initials!
					if visited[v] == 0:						
						nonwater_neighbours.append( (v, u, neighbour_id, id_to_molecule(neighbour_id), s.area()) )

	return visited	


def collect_volumes(time, frame, label, items, some_selection, volumes, to_file):
	# use value key=items[i] to check if volume[i] should be added to 
	# the total value assigned to the key. It is only done if key in some_selection. 
	TMP_RANGE = range(max(some_selection)+1)
	TMP_MIN = min(some_selection)
	TMP_MAX = max(some_selection)
	sum_volumes = [0.0 for d in TMP_RANGE]
	count = [0.0 for d in TMP_RANGE]
	for d, v in zip(items, volumes):
		d = int(d)
		if TMP_MIN <= d <= TMP_MAX:
			if v < 0.0:
				raise Exception("{}: Volume of the cell is NEGATIVE!".format(label))
			sum_volumes[d] += v
			count[d] += 1.0

	data = []
	for d in some_selection:
		if count[d] > 0:
			data.append("{} {} {}".format(sum_volumes[d], count[d], sum_volumes[d] / count[d]))
		else:
			data.append("0 0 0")
	to_file.write("{} {} {}\n" .format(time, frame, " ".join(data)))	
	##########################
	# END of collect_volumes #
	##########################	


def write_file_info(to_file, molecules_keys, layers):
	# output info on the file
	to_file.write("# this file was generated on {} with the following command:\n".format(jobstart))
	to_file.write("# \n")
	to_file.write("#     {}\n".format(" ".join(sys.argv)))
	to_file.write("# \n")
	# output legend
	to_file.write("# LEGEND: \n")
	to_file.write("# time frame              obvious\n")
	to_file.write("# volume_layer_N          total volume in a given layer (or where sum of layers is equal to a given value)\n")	
	to_file.write("# count_layer_N           molecule count in ... (as above)\n")	
	to_file.write("# avgvol_layer_N          per molecule, as there might be different number of molecules in each layer and we want to compare, in ... (as above)\n")
	to_file.write("# \n")
	to_file.write("# This file was generated for job '{}' and molecule(s): '{}'\n".format(jobname, ", ".join(molecules_keys)))
	to_file.write("# \n")
	# output headers
	to_file.write("# (gnuplot column numbers for help)\n")
	to_file.write("# 1    2     ")	
	for i in range(3*len(layers)):
		to_file.write("{}".format(i+3) + (" " * (len("               ") - len("{}".format(i+3)))))
	to_file.write("\n")	
	to_file.write("# time frame ")
	for i in layers:
		to_file.write("volume_layer_{} ".format(i))	
		to_file.write("count_layer_{}  ".format(i))	
		to_file.write("avgvol_layer_{} ".format(i))
	to_file.write("\n")	

# contains the cross water layer with sums of distances (BFS algo) equal something. 
ccw_analys_file = file("{0}/common_crosswater_layers_analysis.dat".format(outf_basename, jobname), "w")
write_file_info(ccw_analys_file, ["(ALL):"] + [molecule_key for molecule_key in molecules], COMMON_CROSS_LAYERS)
#ccw_tjson_file = file("{0}/{1}_common_crosswater.tjson".format(outf_basename, jobname), "w")
common_cross_water = {
	'analys_file': ccw_analys_file, 
#	'detail_file': ccw_tjson_file, 
}
# generate as many analysis files as there are different kinds of molecules in the setup
for molecule_key in molecules:	
	# write here analysis data per frame
	analys_file = file("{0}/{2}_hydbond_layers_analysis.dat".format(outf_basename, jobname, molecule_key), "w")
	write_file_info(analys_file, [molecule_key], COLLECT_LAYERS)
	molecules[molecule_key]['analys_file'] = analys_file

	for d in COLLECT_LAYERS:
		layers_file = file("{0}/{2}_hydbond_layer_L{3}.tjson".format(outf_basename, jobname, molecule_key, d), "w")
		molecules[molecule_key]['layers_L{}_file'.format(d)] = layers_file	

	# we will do cross water layers with other molecules
	molecules[molecule_key]['cross_water'] = {}
	molecules[molecule_key]['areas'] = {}
	for other_key in molecules:	
		# BEWARE: taka sama kolejnosc musi byc jak nizej przy liczeniu
		# BAD CODING :( powinien byc wspolny predykat.
		if other_key > molecule_key:
			cross_file = file("{0}/{2}_crosswater_{3}_analysis.dat".format(outf_basename, jobname, molecule_key, other_key), "w")
			write_file_info(cross_file, [molecule_key, other_key], CROSS_LAYERS)
#			cross_tjson_file = file("{0}/{1}_{2}_crosswater_{3}.tjson".format(outf_basename, jobname, molecule_key, other_key), "w")

			molecules[molecule_key]['cross_water'][other_key] = { 
				'analys_file': cross_file, 
#				'detail_file': cross_tjson_file,
			}

			areas_file = file("{0}/{2}_areas_{3}_analysis.dat".format(outf_basename, jobname, molecule_key, other_key), "w")
			# output info on the file
			areas_file.write("# this file was generated on {} with the following command:\n".format(jobstart))
			areas_file.write("# \n")
			areas_file.write("#     {}\n".format(" ".join(sys.argv)))
			areas_file.write("# \n")
			# output legend
			areas_file.write("# LEGEND: \n")
			areas_file.write("# time frame              obvious\n")
			areas_file.write("# area                    common area betwen given molecules (the areas are disjoint by definition, so you can plot total area of contact as a sum)\n")			
			areas_file.write("# min_distance            minimal distance, if defined by the common faces of voronoi, -1 otherwise\n")			
			areas_file.write("# \n")
			areas_file.write("# This file was generated for job '{}' and molecule(s): '{}'\n".format(jobname, ", ".join([molecule_key, other_key])))
			areas_file.write("# \n")
			# output headers
			areas_file.write("# (gnuplot column numbers for help)\n")
			areas_file.write("# 1    2     3    4 \n")
			areas_file.write("# time frame area min_distance")
			areas_file.write("\n")				
			molecules[molecule_key]['areas'][other_key] = { 'analys_file': areas_file, }			

volume_analys_file = file("{0}/molecule_volume_analysis.dat".format(outf_basename, jobname), "w")
volume_analys_file.write("# this file was generated on {} with the following command:\n".format(jobstart))
volume_analys_file.write("# \n")
volume_analys_file.write("#     {}\n".format(" ".join(sys.argv)))
volume_analys_file.write("# \n")
# output legend
volume_analys_file.write("# LEGEND: \n")
volume_analys_file.write("# time frame              obvious\n")
volume_analys_file.write("# volume_[mol]            volume of the given molecule\n")			
volume_analys_file.write("# exarea_[mol]            area of the contact with water of the molecule\n")			
volume_analys_file.write("# inarea_[mol]            area of the contact with other molecules\n")			
volume_analys_file.write("# \n")
volume_analys_file.write("# This file was generated for job '{}' and molecule(s): '{}'\n".format(jobname, ", ".join([molecule_key, other_key])))
volume_analys_file.write("# \n")
# output headers
volume_analys_file.write("# (gnuplot column numbers for help)\n")
gphelpers = " ".join([ (str(3 + 3*i+j) + (" "*(7 + len(key))))[0:7+len(key)] for (i, key) in enumerate(molecules) for j in range(3) ])
info_headers = " ".join( ["{}_{}".format(valtype, key) for key in molecules for valtype in ["volume", "exarea", "inarea"] ] )
volume_analys_file.write("# 1    2     {}\n".format(gphelpers))
volume_analys_file.write("# time frame {}\n".format(info_headers))	
molecules_volume = { 'analys_file': volume_analys_file, }
pymol_file = {}
for key in molecules:
	pymol_file[key] = []
	for d in set([0,] + list(COLLECT_LAYERS)):
		pymol_file[key].append(file("{0}/layer_{2}_of_{3}.pdb".format(outf_basename, jobname, d, key), "w"))


##################################################################################
##################################################################################
##################################################################################
# THIS IS OLD SCRIPT THAT PARSES PDB FILE AS A TRAJCTORY #########################
# this is not so optimal, but I must for now stick to it #########################
##################################################################################
##################################################################################
##################################################################################
in_file = file(settings.PDB_TRAJ_FILEPATH, "r")  # read from here
try:
	frame = 0
	time = 0
	hao_conversion = [None, ] # first element is just whatever, becouse we start numbering (id atoms) from 1
	hao_ready = False         # I will use this in future to optimize			
	while True: # iterate over all frames	
		header_lines = []
		footer_lines = []
		solvent_lines = []
		hao_lines = []
		molecules_lines = []
		CURR_ATOM_NO = 0 		  # In this frame
		water_positions = {}
		while True: # read one frame			
			line = in_file.readline()        
			if not line:
				raise StopIteration

			if line.startswith("TITLE"):
				header_lines.append(line)			
				# extract time from title
				tpos = line.find("t=");
				if tpos != -1:
					try:
						time = float(line[tpos+2:].strip())
					except Exception as e:
						time = float(line[tpos+2:].strip().split()[0])
				else:
					# this is unfortunate, but I assume file is consistend between frames 
					# (i.e. they all have t or all of them do not have t)
					time = float(frame)

			dat = pdb.tokenize(line)
			if dat is None: 			
				continue;
			item = pdb.PDBRecord(dat)
			if item.type() == "CRYST1":	# create new geometry	
				header_lines.append(line)			
				(a, b, c) = (item["a"], item["b"], item["c"])
				print "CRYST1: Creating geometry and container... ", a, b, c, 
				geometry = OrthogonalGeometry(item["a"], item["b"], item["c"], True, True, True)				
				container = Container(geometry)				
				# this_frame_box = [0.5*a, 0.5*b, 0.5*c]  # TODO: why this one?
				this_pbc_box = [a, b, c]  # pbc_box requires full lengths!
				print "OK"
			elif item.type() == "ATOM":					
				# some bookkeeping for drawing purposes
				if settings.is_solvent(item):
					solvent_lines.append(line)
				else:
					molecules_lines.append(line)					

				if is_hydrogen(item):
					# we allways add hydrogen atoms to preceeding HAO atoms for water (we assume this organization of PDB)
					# for other molecules we do not have this knowledge... (it is more complicated in PDB)
					(atom_id, x, y, z, r) = item.as_particle()
					water_positions[atom_id] = (x, y, z)					
					if settings.is_water(item):
						hao_lines[-1].append(line)						
					# add hydrogens only if HAO (Heavy Atoms Only) is set to false...	
					if (not settings.is_solvent(item) and not mols_hao) or (settings.is_solvent(item) and not solv_hao):
						container.add(atom_id, x, y, z, HYDROGEN_RADIUS) # for now fix r (hydrogen)
				else:
					# all other particles goes into Voronoi diagram
					(atom_id, x, y, z, r) = item.as_particle()
					if settings.is_water(item):
						water_positions[atom_id] = (x, y, z)					
					container.add(atom_id, x, y, z, OTHER_RADIUS) # for now fix r (carbon and others)					
					if not hao_ready:
						print "adding", atom_id, x, y, z, r
						while len(hao_conversion) <= atom_id:
							hao_conversion.append(None)
						hao_conversion[atom_id] = CURR_ATOM_NO
					CURR_ATOM_NO += 1
					# we allways append HAO atoms as lists - we will include Hydrogens later
					# this is ugly, but... 
					hao_lines.append([line, ])				
			elif item.type() == "TER":		
				footer_lines.append(line)					
				print "TER from line", line				
			elif item.type() == "ENDMDL":
				footer_lines.append(line)					
				print "ENDMDL from line", line
				break
			else:
				header_lines.append(line)			
				print "other item type", item.type(), "from line", line
		# OK, NOW WE HAVE PARSED NEXT PDB FRAME, we can do computations.

		hao_ready = True   # After first iteration we have hao_conversion ready for all frames.
		# print hao_conversion
		frame += 1
		if frame > settings.MAX_ANALYSIS_FRAMES:
			break

		###############################
		### HERE WE DO COMPUTATIONS ###
		###############################
		NUM_CELLS = container.size()
		diagram = Diagram(container)
		volumes = [-1.0 for i in range(NUM_CELLS)]				
		for molecule, moldata in molecules.iteritems():
			nonwater_neighbours = []
			water_neighbours = []
			moldata["volume"] = 0.0
			moldata["exarea"] = 0.0
			moldata["inarea"] = 0.0
			moldata["distances"] = {}
			moldata["common_area"] = {}
			initials = [hao_conversion[iid] for iid in range(moldata['def']['lo'], moldata['def']['up']+1) if hao_conversion[iid] is not None]			
			depths = run_BFS(NUM_CELLS, diagram, initials, max(COLLECT_LAYERS), volumes, nonwater_neighbours, water_neighbours)

			# output graphics for pymol
			# if pymol_file[molecule]:
			# 	# we write only once (one frame, as PyMol cannot handle dynamically changing collections)
			# 	for d in set([0,] + list(COLLECT_LAYERS)):
			# 		pymol_file[molecule][d].write("".join(header_lines))
			# 	pymol_file[molecule][0].write("".join(molecules_lines))
			# 	for d, item_num in zip(depths, range(container.size())):
			# 		if d in COLLECT_LAYERS:
			# 			pymol_file[molecule][d].write("".join(hao_lines[item_num]))
			# 	for d in set([0,] + list(COLLECT_LAYERS)):
			# 		pymol_file[molecule][d].write("".join(footer_lines))
			# 		pymol_file[molecule][d].close()
			# 	pymol_file[molecule] = None

			for d in set([0,] + list(COLLECT_LAYERS)):
				pymol_dirpath = "{0}/{1}_L{2}".format(LAYERS_DIRPATH, molecule, d)
				pymol_filepath = "{0}/{4}_L{3}_f{1}_t{2}.pdb".format(pymol_dirpath, frame, time, d, molecule)
				mkdir_p(pymol_dirpath)
				pymol_file = file(pymol_filepath, "w")
				pymol_file.write("".join(header_lines))
				for td, item_num in zip(depths, range(container.size())):
					if td == d:
						pymol_file.write("".join(hao_lines[item_num]))
				pymol_file.write("".join(footer_lines))
				pymol_file.close()

			water_layers = { d: list() for d in COLLECT_LAYERS }
			for num, d in enumerate(depths):
				if d in COLLECT_LAYERS:
					item = container.get(num)
					if settings.is_water(item):
						# append item id (atom ID) (see that only O atoms are appended)
						water_layers[d].append(item[0])

			for d in COLLECT_LAYERS:
				of = molecules[molecule]['layers_L{}_file'.format(d)]
				waters = [(oid, water_positions[oid], water_positions[oid+1], water_positions[oid+2]) for oid in water_layers[d]]
				of.write("{} {}\n".format(time, json.dumps(waters)))

			moldata['depths'] = depths
			collect_volumes(time, frame, "COLLECT_LAYERS", depths, COLLECT_LAYERS, volumes, moldata["analys_file"])		
			common_area = 0.0
			for vno, uno, neighbour_id, neighbour_molecule, side_area in nonwater_neighbours:
				if neighbour_molecule != molecule:					
					if neighbour_molecule not in moldata["distances"]:
						moldata["distances"][neighbour_molecule] = []
					if neighbour_molecule not in moldata["common_area"]:
						moldata["common_area"][neighbour_molecule] = 0.0						
					di = closest_distance(to_coords(container.get(vno)), to_coords(container.get(uno)), this_pbc_box)
					moldata["distances"][neighbour_molecule].append(di)
					moldata["common_area"][neighbour_molecule] += side_area
					moldata["inarea"] += side_area

			for vno, uno, neighbour_id, neighbour_molecule, side_area in water_neighbours:
				moldata["exarea"] += side_area

			for d,v in zip(depths, volumes):
				if d == 0:
					if v < 0.0:
						raise Exception("In molecule volume, volume is negative!")
					moldata["volume"] += v

		# collect common things for all molecules, after all computations and collections done.		
		molecules_volume['analys_file'].write("{} {}".format(time, frame))
		for this_molecule, this_moldata in molecules.iteritems():			
			for other_molecule, other_moldata in molecules.iteritems():
				# correct the 0-th depth - the molecule itself, so that 0 will not be used in collecting water crosslayers				
				# THIS IS A FAST HOT FIX! NOT SO NICE
				for i, d in enumerate(this_moldata['depths']):
					if d == 0:
						this_moldata['depths'][i] = 2 * MINUS_INFTY
				for i, d in enumerate(other_moldata['depths']):
					if d == 0:
						other_moldata['depths'][i] = 2 * MINUS_INFTY
				# BEWARE: taka sama kolejnosc musi byc jak wyzej przy definicji plikow wyjsciowych...
				# BAD CODING :( powinien byc wspolny predykat.
				if other_molecule > this_molecule:
					f = this_moldata['cross_water'][other_molecule]['analys_file']
					items = map(sum, zip(this_moldata['depths'], other_moldata['depths']))
					collect_volumes(time, frame, "CROSS_LAYERS", items, CROSS_LAYERS, volumes, f)

					fa = this_moldata['areas'][other_molecule]['analys_file']
					common_area = this_moldata['common_area'].get(other_molecule, 0.0)
					min_distance = min(this_moldata['distances'].get(other_molecule, [-1.0,]))
					fa.write("{} {} {} {}\n".format(time, frame, common_area, min_distance))

					for d in CROSS_LAYERS:
						pymol_dirpath = "{0}/{1}_crosslayer_{2}_L{3}".format(LAYERS_DIRPATH, this_molecule, other_molecule, d)
						pymol_filepath = "{0}/{4}cross{5}_L{3}_f{1}_t{2}.pdb".format(pymol_dirpath, frame, time, d, this_molecule, other_molecule)
						mkdir_p(pymol_dirpath)
						pymol_file = file(pymol_filepath, "w")
						pymol_file.write("".join(header_lines))
						for td, item_num in zip(items, range(container.size())):
							if td == d:
								pymol_file.write("".join(hao_lines[item_num]))
						pymol_file.write("".join(footer_lines))
						pymol_file.close()

			molecules_volume['analys_file'].write(" {} {} {}".format(this_moldata["volume"], this_moldata["exarea"], this_moldata["inarea"]))
		molecules_volume['analys_file'].write("\n")

		items = map(sum, zip(*[moldata['depths'] for molecule, moldata in molecules.iteritems()]))
		collect_volumes(time, frame, "COMMON_CROSS_LAYERS", items, COMMON_CROSS_LAYERS, volumes, common_cross_water['analys_file'])

		##############
		# NEXT FRAME #
		##############

except StopIteration as e:

	for molecule, moldata in molecules.iteritems():
		for d in COLLECT_LAYERS:
			molecules[molecule]['layers_L{}_file'.format(d)].close()

	print "Finished"
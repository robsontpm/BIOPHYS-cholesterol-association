# coding: utf-8
import sys
sys.path.append('.')
import os

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

# IMPORT SETTINGS
import importlib
import datetime
if len(sys.argv) < 2:
	print "Usage: python orient.py [path.to.conf] {backup_mode}"
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
	if backup_mode not in ["backup", "override"]:
		backup_mode = "backup"
# END IMPORT SETTINGS.

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


def veclen2(x):
	return sum(map(lambda y: float(y)**2, x))


def veclen(x): 
	return sqrt(veclen2(x))	


def normalize(x):
	d = veclen(x)
	if d > 0:
		return map(lambda y: float(y)/d, x)
	else:
		raise Exception("Could not normalize (0,0,0)!")


def crossproduct(x, y):
	return [
		x[1]*y[2]-x[2]*y[1],
		x[2]*y[0]-x[0]*y[2],
		x[0]*y[1]-x[1]*y[0],
	]


def dotproduct(x, y):
	return sum(map(lambda z: z[0]*z[1], zip(x, y)))


def plusv(x, c, y):
	# returns x + c*y, c scalar
	return [
		x[0] + c * y[0],
		x[1] + c * y[1],
		x[2] + c * y[2],
	]	


def minusv(x, c, y):	
	# returns x - c*y, c scalar
	return plusv(x, -c, y)


def dist2(x, y):
	return sum(map(lambda z: (z[0]-z[1])**2, zip(x, y)))


def dist(x, y):
	return sqrt(x, y)


# FOR CHol we have
# def_atom1 = pos of C6
# def_atom2 = pos of C10
# def_atom3 = pos of C13
# on output we have:
# C10-C13 direction is preserved (but normalized)
# C6-C10 direction is ortonormalized to C10-C13
def compute_coords(def_atom1, def_atom2, def_atom3):
	# compute vectors normalized to (0,0,0) in def_atom2
	v = [
		def_atom1[0] - def_atom2[0],
		def_atom1[1] - def_atom2[1],
		def_atom1[2] - def_atom2[2],
	]
	u = [
		def_atom3[0] - def_atom2[0],
		def_atom3[1] - def_atom2[1],
		def_atom3[2] - def_atom2[2],
	]	
	u = normalize(u)
	v = minusv(v, dotproduct(u, v), u)
	v = normalize(v)
	# vectors v, u defines a plane, compute w - third orthonormal vector
	w = crossproduct(v, u)
	return u, v, w


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


def sinsign(value):
	return 1.0 * (value < 0)


# generate as many analysis files as there are different kinds of molecules in the setup
for molecule_key in molecules:	
	molecules[molecule_key]['angles'] = {}
	for other_key in molecules:	
		# BEWARE: taka sama kolejnosc musi byc jak nizej przy liczeniu
		# BAD CODING :( powinien byc wspolny predykat.
		if other_key > molecule_key:
			angles_file = file("{0}/{1}_{2}_angles_{3}_analysis.dat".format(outf_basename, jobname, molecule_key, other_key), "w")
			# output info on the file
			angles_file.write("# this file was generated on {} with the following command:\n".format(jobstart))
			angles_file.write("# \n")
			angles_file.write("#     {}\n".format(" ".join(sys.argv)))
			angles_file.write("# \n")
			# output legend
			angles_file.write("# LEGEND: \n")
			angles_file.write("# time frame              obvious\n")
			angles_file.write("# u v w                   angle (radians) of the defining vector (w is just normal to the defining plane) between two molecules. (0,0,0) is a perfect align\n")
			angles_file.write("#                         u, v are given by the defining atoms, u: 2->1, v: 2->3, in the computation, v is preserved, u is orthogonalised to v,\n")
			angles_file.write("# ua va                   angle (radians) of the defining vector, as full range 0-2pi (decided based on the side of the definining plane u,v), \n")
			angles_file.write("# d                       closest distance from the three given later\n")
			angles_file.write("# d0 d1 d2                distance to the closest images between coresponding defining atoms, the lower three numbers the better align\n")
			angles_file.write("# \n")
			angles_file.write("# This file was generated for job '{}' and molecule(s): '{}'\n".format(jobname, ", ".join([molecule_key, other_key])))
			angles_file.write("# \n")
			# output headers
			angles_file.write("# (gnuplot column numbers for help)\n")
			angles_file.write("# 1    2     3 4 5 6  7  8 9  10 11 \n")
			angles_file.write("# time frame u v w ua va d d0 d1 d2")
			angles_file.write("\n")				
			molecules[molecule_key]['angles'][other_key] = { 'analys_file': angles_file, }			


def_out_traj_lo = min(item['def']['lo'] for key, item in molecules.iteritems())
def_out_traj_up = max(item['def']['up'] for key, item in molecules.iteritems())
out_traj_filepath = os.path.join(settings.OUTPUT_DIRPATH, "orient_traj.pdb")
out_traj_file = file(out_traj_filepath, "w")

out_mols_traj_filepath = os.path.join(settings.OUTPUT_DIRPATH, "orient_traj_mols.pdb")
out_mols_traj_file = file(out_mols_traj_filepath, "w")
out_orients_traj_filepath = os.path.join(settings.OUTPUT_DIRPATH, "orient_traj_orients.pdb")
out_orients_traj_file = file(out_orients_traj_filepath, "w")

#orients_scale = 15.0
orients_scale = 1.0
sel_filepath = os.path.join(settings.OUTPUT_DIRPATH, "orient_traj.pymol")
sel_file = file(sel_filepath, "w")
sel_orients_filepath = os.path.join(settings.OUTPUT_DIRPATH, "orient_traj_orients.pymol")
sel_orients_file = file(sel_orients_filepath, "w")
for imol, item in enumerate(molecules):
	for a, color in enumerate(["white", "red", "green", "blue"]):
		sel_file.write("sel imol{}a{}, id {}\n".format(imol, a, def_out_traj_up + 1 + imol * 4 + a))
		sel_file.write("color {}, imol{}a{}\n".format(color, imol, a))
		sel_orients_file.write("sel imol{}a{}, id {}\n".format(imol, a, 1 + imol * 4 + a))
		sel_orients_file.write("color {}, imol{}a{}\n".format(color, imol, a))		
sel_file.close()
sel_orients_file.close()



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
		CURR_ATOM_NO = 0 		  # In this frame
		container = []
		while True: # read one frame			
			line = in_file.readline()        
			if not line:
				raise StopIteration

			if line.startswith("TITLE"):
				out_traj_file.write(line)	
				out_orients_traj_file.write(line)			
				out_mols_traj_file.write(line)	
				# extract time from title
				tpos = line.find("t=");
				if tpos != -1:
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
				out_traj_file.write(line)
				out_orients_traj_file.write(line)
				out_mols_traj_file.write(line)	
				(a, b, c) = (item["a"], item["b"], item["c"])			
				this_frame_box = [0.5*a, 0.5*b, 0.5*c]  # TODO: why this one?
				this_pbc_box = [a, b, c]  # pbc_box requires full lengths!
			elif item.type() == "ATOM":			
				# add hydrogens only if HAO (Heavy Atoms Only) is set to false...	
				(atom_id, x, y, z, r) = item.as_particle()
				if is_hydrogen(item):
					if (not settings.is_solvent(item) and not mols_hao) or (settings.is_solvent(item) and not solv_hao):						
						container.append( (atom_id, x, y, z, HYDROGEN_RADIUS) ) # for now fix r (hydrogen)
				else:
					# all other particles goes into Voronoi diagram
					container.append( (atom_id, x, y, z, OTHER_RADIUS) ) # for now fix r (carbon and others)					
					if not hao_ready:
						while len(hao_conversion) <= atom_id:
							hao_conversion.append(None)
						hao_conversion[atom_id] = CURR_ATOM_NO
					CURR_ATOM_NO += 1

				if def_out_traj_lo <= atom_id <= def_out_traj_up:
					out_traj_file.write(line)
					out_mols_traj_file.write(line)	
			elif item.type() == "TER":			
				print "TER from line", line			
				break
			elif item.type() == "ENDMDL":
				print "ENDMDL from line", line
			else:
				out_traj_file.write(line)
				out_orients_traj_file.write(line)
				out_mols_traj_file.write(line)	
				print "other item type", item.type(), "from line", line
		# OK, NOW WE HAVE PARSED NEXT PDB FRAME, we can do computations.

		hao_ready = True   # After first iteration we have hao_conversion ready for all frames.
		frame += 1
		if frame > settings.MAX_ANALYSIS_FRAMES:
			break

		for imol, (this_molecule, this_moldata) in enumerate(molecules.iteritems()):
			this_da_ids = this_moldata['def']['plane_atoms_ids']
			this_da = [
				to_coords(container[hao_conversion[this_da_ids[0]]]), 
				to_coords(container[hao_conversion[this_da_ids[1]]]), 
				to_coords(container[hao_conversion[this_da_ids[2]]]),
			]
			this_u, this_v, this_w = compute_coords(this_da[0], this_da[1], this_da[2])

			local_u = plusv(plusv(this_u, orients_scale-1.0, this_u), 1.0, this_da[1])   # wektorowo: orients_scale * this_u + this_da[1]  (skala, zeby cos bylo widac!, mozna ustalic na poczatku skryptu)
			local_v = plusv(plusv(this_v, orients_scale-1.0, this_v), 1.0, this_da[1])
			local_w = plusv(plusv(this_w, orients_scale-1.0, this_w), 1.0, this_da[1])
			out_traj_file.write(pdb.hetatm_line(def_out_traj_up + 1 + imol * 4 + 0, this_da[1][0], this_da[1][1], this_da[1][2], 2.0))
			out_traj_file.write(pdb.hetatm_line(def_out_traj_up + 1 + imol * 4 + 1, local_u[0], local_u[1], local_u[2], 1.0))
			out_traj_file.write(pdb.hetatm_line(def_out_traj_up + 1 + imol * 4 + 2, local_v[0], local_v[1], local_v[2], 1.0))
			out_traj_file.write(pdb.hetatm_line(def_out_traj_up + 1 + imol * 4 + 3, local_w[0], local_w[1], local_w[2], 1.0))
			out_traj_file.write(pdb.conect_line(def_out_traj_up + 1 + imol * 4 + 0, def_out_traj_up + 1 + imol * 4 + 1))
			out_traj_file.write(pdb.conect_line(def_out_traj_up + 1 + imol * 4 + 0, def_out_traj_up + 1 + imol * 4 + 2))
			out_traj_file.write(pdb.conect_line(def_out_traj_up + 1 + imol * 4 + 0, def_out_traj_up + 1 + imol * 4 + 3))

			# TODO: DRY :(
			out_orients_traj_file.write(pdb.hetatm_line(1 + imol * 4 + 0, this_da[1][0], this_da[1][1], this_da[1][2], 2.0))
			out_orients_traj_file.write(pdb.hetatm_line(1 + imol * 4 + 1, local_u[0], local_u[1], local_u[2], 1.0))
			out_orients_traj_file.write(pdb.hetatm_line(1 + imol * 4 + 2, local_v[0], local_v[1], local_v[2], 1.0))
			out_orients_traj_file.write(pdb.hetatm_line(1 + imol * 4 + 3, local_w[0], local_w[1], local_w[2], 1.0))
			out_orients_traj_file.write(pdb.conect_line(1 + imol * 4 + 0, 1 + imol * 4 + 1))
			out_orients_traj_file.write(pdb.conect_line(1 + imol * 4 + 0, 1 + imol * 4 + 2))
			out_orients_traj_file.write(pdb.conect_line(1 + imol * 4 + 0, 1 + imol * 4 + 3))			

			for other_molecule, other_moldata in molecules.iteritems():
				if other_molecule > this_molecule:
					# TODO: FIX: DRY
					other_da_ids = other_moldata['def']['plane_atoms_ids']
					other_da = [
						to_coords(container[hao_conversion[other_da_ids[0]]]), 
						to_coords(container[hao_conversion[other_da_ids[1]]]), 
						to_coords(container[hao_conversion[other_da_ids[2]]]),
					]
					other_u, other_v, other_w = compute_coords(other_da[0], other_da[1], other_da[2])
					cos_v = dotproduct(this_v, other_v)
					cos_u = dotproduct(this_u, other_u)
					cos_w = dotproduct(this_w, other_w)
					sign_sin_v = sinsign(dotproduct(this_w, other_v))
					sign_sin_u = sinsign(dotproduct(this_w, other_u))
					simple_angle_u = acos(cos_u)  #acos(dotproduct(this_v, other_v))
					simple_angle_v = acos(cos_v)  #acos(dotproduct(this_u, other_u))
					simple_angle_w = acos(cos_w)  #acos(dotproduct(this_w, other_w))
					full_angle_u = simple_angle_u + sign_sin_v * pi
					full_angle_v = simple_angle_v + sign_sin_u * pi
					d0 = closest_distance(this_da[0], other_da[0], this_pbc_box)
					d1 = closest_distance(this_da[1], other_da[1], this_pbc_box)
					d2 = closest_distance(this_da[2], other_da[2], this_pbc_box)
					f = this_moldata['angles'][other_molecule]['analys_file']
					f.write("{} {} {} {} {} {} {} {} {} {} {}\n".format(
						time, frame, 
						simple_angle_u, simple_angle_v, simple_angle_w, 
						full_angle_u, full_angle_v,						
						min(d0, d1, d2), d0, d1, d2
					))

		out_traj_file.write("TER\n")	
		out_traj_file.write("ENDMDL\n")	
		out_orients_traj_file.write("TER\n")	
		out_orients_traj_file.write("ENDMDL\n")			
		out_mols_traj_file.write("TER\n")	
		out_mols_traj_file.write("ENDMDL\n")				

		##############
		# NEXT FRAME #
		##############

except StopIteration as e:
	print "Finished"


out_traj_file.close()
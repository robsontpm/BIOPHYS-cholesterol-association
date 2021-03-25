# coding: utf-8
import sys
sys.path.append('.')
import os

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


# # IMPORT SETTINGS
if (len(sys.argv) < 2):
	print "Computes the first hydratation layer statistics."
	print ""
	print "Usage: python voro3d.py [path.to.conf] {no-backup}"
	print ""
	print "  [path.to.conf] is a python module path to a file containing configuration.";
	print "                 see sample files in jobs/* for reference."
	print "  {no-backup}    optional. If theres is no-backup in the end, "
	print "                 then the output directory will not be backuped."
	print "                 By default for safety, we always do backup."
	print "                 No backup is good for recomputing or computing to "
	print "                 common directory by many scripts."
	sys.exit(-1)

job = helper.Job(sys.argv)
outf_basename = job.settings.OUTPUT_DIRPATH

# common job.settings.
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
try: solv_hao = job.settings.SOLVENT_HAO
except: job.log("# WARNING: using default SOLVENT_HAO = {}\n".format(solv_hao))
try: mols_hao = job.settings.MOLECULES_HAO
except: job.log("# WARNING: using default MOLECULES_HAO = {}\n".format(mols_hao))
try: max_analys_frames = job.settings.MAX_ANALYSIS_FRAMES
except: job.log("# WARNING: using default MAX_ANALYSIS_FRAMES = {}\n".format(max_analys_frames))
try: max_voro_frames = job.settings.MAX_VORO_FRAMES 
except: job.log("# WARNING: using default MAX_VORO_FRAMES = {}\n".format(max_voro_frames))
try: molecules = job.settings.MOLECULES; 
except: job.log("# WARNING: using default MOLECULES (CHOL) = {}\n".format(str(molecules)))

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


in_file = file(job.settings.PDB_TRAJ_FILEPATH, "r")  # read from here
vorobasedir = "{0}/voro/".format(outf_basename)
if not os.path.exists(vorobasedir):
	os.mkdir(vorobasedir)

# generate as many analysis files as there are different kinds of molecules in the setup
for molecule_key in molecules:	
	voro_filename_template = "{0}/{1}_{2}_compute3d_voro_".format(vorobasedir, job.jobname, molecule_key) + "{:04d}.pdb";
	molecules[molecule_key]['voro_filename_template'] = voro_filename_template

try:
	frame = 0
	time = 0
	while True: # iterate over all frames	

		# THIS IS AN OLD SCRIPT SO IT NEED TO PARSE PDB BY ITS OWN...
		# PARSING STARTS HERE...
		lines = []
		particles = []
		items = []
		header = []
		bottom = []

		while True:
			line = in_file.readline()        
			if not line:
				raise StopIteration

			if line.startswith("TITLE"):
				tpos = line.find("t=");
				if tpos != -1:
					time = float(line[tpos+2:].strip())
				else:
					time = float(frame)

			dat = pdb.tokenize(line)
			if dat is None: 
				header.append(line)
				continue;
			item = pdb.PDBRecord(dat)
			if item.type() == "CRYST1":	# create new geometry
				print "Creating geometry... "
				(a, b, c) = (item["a"], item["b"], item["c"])
				geometry = OrthogonalGeometry(item["a"], item["b"], item["c"], True, True, True)				
				print "Creating container... "
				container = Container(geometry)				
				print "Geometry and Container OK... "				
				header.append(line)
			elif item.type() == "ATOM":			
				lines.append(line)
				items.append(item)	

				# add hydrogens only if HAO (Heavy Atoms Only) is set to false...	
				if is_hydrogen(item):
					if (not job.settings.is_solvent(item) and not mols_hao) or (job.settings.is_solvent(item) and not solv_hao):
						(atom_id, x, y, z, r) = item.as_particle()
						container.add(atom_id, x, y, z, HYDROGEN_RADIUS) # for now fix r (hydrogen)
				else:
					# all other particles goes into Voronoi diagram
					(atom_id, x, y, z, r) = item.as_particle()
					container.add(atom_id, x, y, z, OTHER_RADIUS) # for now fix r (carbon and others)
			elif item.type() == "TER":				
				print "TER from line", line
				break
			elif item.type() == "ENDMDL":
				bottom.append(line)
			else:
				print "other item type", item.type(), "from line", line
				header.append(line)
		# OK, NOW WE HAVE PARSED NEXT PDB FRAME, we can do computations.

		# tu sie buduje Voronoi
		diagram = Diagram(container)

		# to jest hack na animacje: tworze tyle plikow PDB ile klatek do zapisania voronoia
		# po czym dopisuje puste klatki. Tylko w klatce frame pisze w pliku frame wspolrzedne
		# pozniej w pymolu trzeba otworzyc wszystkie na raz. Jest to chyba spowodowane problemami z 
		# numerowaniem wierzchołkow.
		for key, moldata in molecules.iteritems():

			# this saves a PDBs that can be exported to PyMol as animation
			# we need several empty models, then a current model, then again empty models. 
			for kkk in range(0, max_voro_frames):		
				if frame == 0:
					# rewrite file new (if exists) if frame == 0
					voro_file = file(moldata['voro_filename_template'].format(kkk), "w")
				elif frame < max_voro_frames:
					# for other frames - append to that file
					voro_file = file(moldata['voro_filename_template'].format(kkk), "a")

				# save headers!
				if frame < max_voro_frames:
					for l in header:
						voro_file.write(l);

				# flush data!
				voro_file.close()

			# open file corresponding to current frame to write Voronoi diagram here.
			if frame < max_voro_frames:
				voro_file = file(moldata['voro_filename_template'].format(frame), "a")
				moldata['voro_file'] = voro_file

			# since we are inside a loop over molecules, I can initialize important variables
			# that are local to this frame diagram and individual for each molecule
			moldata['base_no'] = 0

		# tutaj zaczynamy iteracje po diagramie
		for i, cell in enumerate(diagram):
			this = container.get(i)
			molkey = id_to_molecule(this[0])
			if not molkey:			
				break

			# we have identified malecule that this cell belongs to, so we
			# will use its 'local' variables from now on. 
			moldata = molecules[molkey]		
			# voro_vertices might be local, as it is used only for drawing voronoi cell for this atom
			voro_vertices = []					
			for side in cell.sides:
				other = container.get(side.neighbour)
				neighbour_id = other[0]
				if not(moldata['def']['lo'] <= neighbour_id <= moldata['def']['up']): 
					# voro_vertices used to draw Voro later
					voro_vertices += side.as_list()

			# TODO: ten IF dodałem, powinien zadziałac, bo do niczego innego
			# tutaj tego malowania scian nie uzwam....
			if frame < max_voro_frames:
				# plus_no is local for this loop and not used anymore after exit, so it might be local
				plus_no = 0
				for v in cell.vertices:
					plus_no += 1						
					(no, x, y, z) = v.as_tuple()
					# we output only vertices for sides that are connected to water
					# therefore we do not produce 'inside' contact surfaces.
					if no in voro_vertices:
						# TODO: two lines below are if we had align3d together with voro3d.
						# TODO: think - maybe save Transformations from align3d and use them here...
						# if moldata['trans']:
						# 	[x, y, z] = moldata['trans'].apply([x,y,z])
						moldata['voro_file'].write(pdb.hetatm_line(moldata['base_no'] + no, x, y, z, 1.0))
				
				for e in cell.edges:
					if e.first() in voro_vertices and e.second() in voro_vertices:
						moldata['voro_file'].write(pdb.conect_line(moldata['base_no'] + e.first(), moldata['base_no'] + e.second()))				

				moldata['base_no'] += plus_no

			# now we go to the next cell

		# here we have finished iterating over cells in a diagram.
		# zakoncz wypisywanie klatek voronoia do kolejnych plikow
		# patrz dokumentacja wyzej, dlaczego tak robimy.
		if frame < max_voro_frames:
			for key, moldata in molecules.iteritems():
				# close all files corresponding to current frame to flush buffers		
				if moldata['voro_file']:
					moldata['voro_file'].close()
					moldata['voro_file'] = None
				# now end all models in all frames and all models.
				for kkk in range(0, max_voro_frames):				
					voro_file = file(moldata['voro_filename_template'].format(kkk), "a")			
					for l in bottom:
						voro_file.write(l);
					voro_file.write("TER\n");					


		# here all molecules are computed, we can move to the next frame
		frame += 1
		if frame >= job.settings.MAX_VORO_FRAMES:
			break

except Exception as e:
	job.finish("By exception")
	print e
	traceback.print_exc()
else:
	job.finish("Finished normally")



# coding: utf-8
import sys
sys.path.append('.')
import os

from pydmga.io import pdb
from pydmga.io import helper
import traceback

from math import pi
from math import sqrt
from math import cos
from math import sin
from math import acos
from math import asin

import numpy as np

#######################
# HELPER FUNCTIONS    #
#######################

def rotMatrix(u, theta):
	"""
	Generates 3D rotation maatrix by an angle theta with axis given by u
	"""
	if isinstance(u, (list, tuple)):
		(ai, bj, ck) = u
	else:
		(ai, bj, ck) = (u.item((0,0)), u.item((1,0)), u.item((2,0)))
	q0 = cos(theta/2.0)
	q1 = sin(theta/2.0) * ai
	q2 = sin(theta/2.0) * bj
	q3 = sin(theta/2.0) * ck

	Q = [	[q0**2 + q1**2 - q2**2 - q3**2, 	2.0 * (q1*q2 - q0*q3), 			2.0 * (q1*q3 + q0*q2)			],
			[2.0 * (q2*q1 + q0*q3), 			q0**2 - q1**2 + q2**2 - q3**2,	2.0 * (q2*q3 - q0*q1)			],
			[2.0 * (q3*q1 - q0*q2),				2.0 * (q3*q2 + q0*q1),			q0**2 - q1**2 - q2**2 + q3**2	]	]

	Q = np.matrix(Q)
	return Q


def normalized(v):
	"""
	Fail-safe version of normalization of the vector
	"""
	norm=np.linalg.norm(v)
	if norm==0: 
		return v
	return v/norm	


def rotVectorToVector(q, p):
	"""
	returns matrix of rotation from q to p
	so that Q * q \in Span({p})
	"""
	v = normalized(p)
	w = normalized(q)
	s = np.linalg.norm(np.cross(v, w, axis=0))
	c = v.transpose().dot(w)
	# theta = asin(s)
	theta = acos(c)		
	t = normalized(np.cross(v, w, axis=0))
	Q1 = rotMatrix(t,  theta)
	Q2 = rotMatrix(t, -theta)
	testw1 = Q1*q
	testw2 = Q2*q
	norm1 = np.linalg.norm(testw1 - float((testw1.transpose()).dot(v)) * v)
	norm2 = np.linalg.norm(testw2 - float((testw2.transpose()).dot(v)) * v)
	if (norm1 < norm2):
		return Q1
	else:
		return Q2


class MoleculeTransform:
	"""
	A transformation from a given data frame to the reference frame
	"""
	def __init__(self, pattern, data):
		"""
		pattern is the original coordinate frame 
		data is the current coordinate frame we want to rotate to fit pattern
		coordinate frames are given as a triangle of points - they define a plane in 3D
		then a third vector is defined as perpendicular to this plane. 
		"""
		p1, p2, p3 = [np.array([p]).transpose() for p in pattern]
		q1, q2, q3 = [np.array([p]).transpose() for p in data]
		pat1, pat2, pat3 = p1, p2, p3
		dat1, dat2, dat3 = q1, q2, q3
		self.tp = p3 # initial translation to 0.0 is -tp in reference frame
		self.tq = q3 # initial transaltion to 0.0 is -tq in data set		
		p1 = p1 - p3
		p2 = p2 - p3
		p3 = p3 - p3
		q1 = q1 - q3
		q2 = q2 - q3
		q3 = q3 - q3
		# orthonormalize base (q1, q2) and (p1, p2)
		q2 = normalized(q2)
		p2 = normalized(p2)
		q1 = normalized(q1 - float( (q1.transpose()).dot(q2) ) * q2)
		p1 = normalized(p1 - float( (p1.transpose()).dot(p2) ) * p2)

		# make rotation matrix that takes q2 to Span({p2})
		self.Q = rotVectorToVector(q2, p2)
		# rotate base points
		q2 = self.Q * q2
		q1 = self.Q * q1
		# make rotation matrix that takes Q * q1 to Span({p1})
		self.W = rotVectorToVector(q1, p1)
		self.M = self.W * self.Q

		# now apply in this order: 
		# a) point -= self.tq
		# b) point = self.Q * point
		# b) point = self.W * point
		# b) point += self.tp
		# and we should have two frames aligned

		# test, dots below should be small (around 0)
		tdat1 = self.apply(dat1) - self.tp
		tdat2 = self.apply(dat2) - self.tp
		tdat3 = self.apply(dat3) - self.tp		

		t = normalized(np.cross(p1, p2, axis=0)).transpose()

	def transform(self, point):
		"""
		Apply transform to any other point
		Point must be np.array column vector. 
		"""
		return self.M * (point - self.tq) + self.tp

	def apply(self, point):
		"""
		Apply transform to any other point
		Point can be np.array column vector
		or any python iterable with 3 coordinates.
		Result returned will be in the same format
		"""
		if isinstance(point, (list, tuple)):
			nppoint = np.array([point]).transpose()
			nppoint = self.transform(nppoint)
			nppoint.transpose()
			return [float(p) for p in nppoint]
		else:
			return self.transform(point)


###############################################################
# REAL SCRIPT STARTS HERE #####################################
###############################################################

if __name__ == "__main__":

	# # IMPORT SETTINGS
	if (len(sys.argv) < 2):
		print "Align molecules between frames."
		print ""
		print "Usage: python align3d.py [path.to.conf] {no-backup}"
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

	# this is some initial setting. SHould be overwritten in settings.
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

	# make a dict that translates ATOM_ID to MOLECULE key. 
	id_to_mol = dict()
	for key, moldata in molecules.iteritems():
		for i in range(moldata['def']['lo'], moldata['def']['up']+1):
			id_to_mol[i] = key


	# we define here local helper function that uses static data from settings
	# TODO: this is not the most clean solution...
	def id_to_molecule(i):
		if i in id_to_mol:
			return id_to_mol[i]
		return None


	# we define here local helper function that uses static data from settings
	# TODO: this is not the most clean solution...
	def num_to_molecule(i):
		if (i+1) in id_to_mol:
			return id_to_mol[i+1]
		return None


	in_file = file(job.settings.PDB_TRAJ_FILEPATH, "r")  # read from here
	# generate as many analysis files as there are different kinds of molecules in the setup
	for molecule_key in molecules:	
		 # write here aligned cholesterols
		transf_out_file = job.user_file("{jobname}_{molecule_key}_align3d_trsf.pdb", {"molecule_key": molecule_key})
		molecules[molecule_key]['transf_out_file'] = transf_out_file

		molecules[molecule_key]['prev_coords'] = None
		molecules[molecule_key]['curr_coords'] = None

	try:
		frame = 0
		time = 0
		while True: # iterate over all frames	

			# THIS IS AN OLD SCRIPT SO IT NEED TO PARSE PDB BY ITS OWN...
			# PARSING STARTS HERE...
			lines = []
			particles = []
			items = []
			bottom = []
			
			for key in molecules:
				# if molecules[key]['prev_coords'] are defined then do nothing, else push there curr_cords
				molecules[key]['prev_coords'] = molecules[key]['prev_coords'] or molecules[key]['curr_coords']
				# make current cords empty at each frame
				molecules[key]['curr_coords'] = list()

			while True:
				line = in_file.readline()  
				if not line:
					raise StopIteration

				dat = pdb.tokenize(line)
				if dat is None: 
					# this is unknown type of token
					map(lambda key: molecules[key]['transf_out_file'].write(line), molecules); # write to all...

				item = pdb.PDBRecord(dat)
				if item.type() == "TITLE":
					time = item["time"] or float(frame)
				if item.type() == "CRYST1":	# create new geometry			
					map(lambda key: molecules[key]['transf_out_file'].write(line), molecules); # write to all...
				elif item.type() == "ATOM":			
					# append to current coords if atom is a defining atom for some molecule
					for key in molecules:
						if item["id"] in molecules[key]['def']['plane_atoms_ids']:
							molecules[key]['curr_coords'].append(item.as_coords())

					lines.append(line)
					items.append(item)	
				elif item.type() == "TER":				
					print "Reading ", frame, " model finished!"
					# map(lambda key: molecules[key]['transf_out_file'].write(line), molecules); # write to all...
					# do not write TER becouse we need to write coordinates first!
					bottom.append(line)
					break
				elif item.type() == "ENDMDL":
					# map(lambda key: molecules[key]['transf_out_file'].write(line), molecules); # write to all...
					# do not write becouse we need to write coordinates first!
					bottom.append(line)
				else:
					print "WARNING: Unknown item type", item.type(), "from line", line
					map(lambda key: molecules[key]['transf_out_file'].write(line), molecules); # write to all...
			# OK, NOW WE HAVE PARSED NEXT PDB FRAME, we can do computations.

			# THIS COMPUTES TRANSFORMATION TO A COMMON FRAME (TO SEE HOW THE MOLECULE 'TWISTS')
			# WE WRITE TRANSFORMED MOLECULES AS CONSECUTIVE FRAMES TO A PDB FILE.

			# first, we make transformation to be applied		
			for key, moldata in molecules.iteritems():
				# we need to store moldata['trans'], as it will be used for 
				# transforming voronoi cells later.
				if moldata['prev_coords'] is not None:
					moldata['trans'] = MoleculeTransform(moldata['prev_coords'], moldata['curr_coords'])
				else:
					moldata['trans'] = None

				for i, p in enumerate(items):
					# i will be used to got to the lines...
					if moldata['def']['lo'] <= p["id"] <= moldata['def']['up']:
						if moldata['trans']:
							pp = moldata['trans'].apply(p.as_coords())
						else:
							pp = p.as_coords()
						linea = lines[i][:30]
						lineb = lines[i][54:]
						moldata['transf_out_file'].write(
							"{0}{1}{2}{3}{4}".format(
								linea, 
								pdb.coord_str(pp[0]), 
								pdb.coord_str(pp[1]), 
								pdb.coord_str(pp[2]),
								lineb
							)
						)
				moldata['transf_out_file'].write("\n".join(bottom));
			# END OF COMPUTING AND WRITING TRANSFORMATION

			# here all molecules are computed, we can move to the next frame
			frame += 1

	except Exception as e:
		job.finish("Finished with exception: {}".format(str(e)))
		# for key, moldata in molecules.iteritems():
		# 	if moldata['transf_out_file']:
		# 		moldata['transf_out_file'].close()		
		print "Finished with exception: ", e
		traceback.print_exc()
	else:
		print "Finished without exception"
		job.finish("Finished without exception")

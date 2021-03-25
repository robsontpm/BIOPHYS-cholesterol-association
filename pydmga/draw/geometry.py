"""
This file defines geometry used to read and write data about
balls in 3D, their voronoi diagrams and motion (in time).

it is used by:
 * voro_render.py
 * voro_render_gui.py
 * voro_gen.py
"""
import random
from math import pi, sin, cos, sqrt
import math
import os #os.path.*
import re #regular expressions (split)
import sys #sys.argv

class Position:
	def __init__(self, x = 0.0, y = 0.0, z = 0.0):
		self.x = x
		self.y = y
		self.z = z			

class Ball:
	def __init__(self, id, x = 0.0, y = 0.0, z = 0.0, r = 1.0):
		self.id = id
		#TODO: other data from atom(mol, atom, etc.)?
		self.x = x
		self.y = y
		self.z = z
		self.r = r

class OneFrame:
	"""OneFrame - class of one Frame of animation of voronoi diagram
	
	For now its very simple
	"""
	def __init__(self):
		self.balls = []
		self.voro = []

class Input:
	"""base class for moving voronoi, it allows for navigating through
	the animation. It's like interface - for me to not forgot what methods are needed ;)
	"""
	def first_frame(self):
		raise Exception("Not implemented yet")
		#do nothing
		
	def last_frame(self):
		raise Exception("Not implemented yet")
		#do nothing		
		
	def current_frame(self):
		raise Exception("Not implemented yet")
		#do nothing
	
	def __init__(self):
		self.balls = []
		self.current = OneFrame()
		
	def advance(self):
		raise Exception("Not implemented yet")
		#do nothing
		
	def go(self, frame):
		raise Exception("Not implemented yet")
		#do nothing
		
	def leap(self, offset):
		raise Exception("Not implemented yet")
		#do nothing	
		
class ArrayInput(Input):
	def first_frame(self):
		return 0;
		
	def last_frame(self):
		return 0;		
		
	def current_frame(self):
		return 0;
	
	def __init__(self, balls, cells):
		self.balls = balls
		self.current = OneFrame()
		self.current.balls = balls
		self.current.voro = cells
		
	def advance(self):
		return False
		
	def go(self, frame):
		return False
		
	def leap(self, offset):
		return False
		
		
class MultiFileInput(Input):
	""" Its kind of input that uses multiple files to store information:
	one file with full description
	(in this way its wery simple to program
	"""
	def _make_current_path(self):
		self._coord_file = os.path.join(self._directory, self._basename, self._basename + ('%08d' % self._current_frame) + '.coord' )
		self._voro_file = os.path.join(self._directory, self._basename, self._basename + ('%08d' % self._current_frame) + '.voro' )
		
	def _load_balls(self):
		self.balls = []
		f = open(self._description_file)
		string_box = re.split('\s+', f.readline().strip())
		self.box_x = float(string_box[0])
		self.box_y = float(string_box[1])
		self.box_z = float(string_box[2])
		self.box_periodic_x = int(string_box[3])
		self.box_periodic_y = int(string_box[4])
		self.box_periodic_z = int(string_box[5])
		#TODO: store box
		for line in f:
			#strip \n, split by whitespaces, convert to types...
			string_data = re.split('\s+', line.strip())
			id = int(string_data[0])
			name = string_data[1]
			moleculeId = int(string_data[2])
			symbol = string_data[3]
			x = float(string_data[4]) - self.box_x/2.0
			y = float(string_data[5]) - self.box_y/2.0
			z = float(string_data[6]) - self.box_z/2.0
			r = float(string_data[7])
			#add new ball , TODO: other data?
			self.balls.append(Ball(id, x, y, z, r))
		f.close()
		
	def _load_frame(self):
		new_frame = OneFrame()
		#read new ball positions
		f_coord = open(self._coord_file)
		for line in f_coord:
			#strip \n, split by whitespaces, convert to types...
			string_data = re.split('\s+', line.strip())
			x = float(string_data[0]) - self.box_x/2.0
			y = float(string_data[1]) - self.box_y/2.0
			z = float(string_data[2]) - self.box_z/2.0
			new_frame.balls.append(Position(x, y, z))		
		f_coord.close()
		#read voro structure
		f_voro = open(self._voro_file)
		line = f_voro.readline();
		while line != '':
			num_points = int(line.strip());
			points = []
			for i in range(num_points):
				line = f_voro.readline();
				string_data = re.split('\s+', line.strip())
				x = float(string_data[0]) - self.box_x/2.0
				y = float(string_data[1]) - self.box_y/2.0
				z = float(string_data[2]) - self.box_z/2.0
				if (len(string_data) > 3):
					point_data = eval(string_data[3])					
				else:
					point_data = {'alpha': 0.0}
				points.append( (x,y,z,point_data) );
			line = f_voro.readline();
			num_faces = int(line.strip());
			faces = []
			for i in range(num_faces):
				line = f_voro.readline();				
				string_data = re.split('\s+', line.strip())
				if (string_data[0].strip('{}') != string_data[0]):
					poly_data = eval(string_data[0])
					string_data.pop(0)
				else:
					poly_data = {'alpha': 0.0}
								
				i = 0
				n = len(string_data)
				poly = []
				while (i < n):
					if (string_data[i].strip('{}') != string_data[i]):
						poly[-1] = (poly[-1][0], eval(string_data[i]))
					else:
						poly.append( (int(string_data[i]), {'alpha': 0.0}) )
					i = i + 1
				
				faces.append((poly,poly_data))
				
			new_frame.voro.append( (points, faces) )
			line = f_voro.readline()
		f_voro.close()
		#store it
		self.current = new_frame

	def __init__(self, repository_path, start_frame = 0, end_frame = None):
		self._repository_path = repository_path
		
		(self._directory, self._basename) = os.path.split(repository_path)
		if self._basename == '':
			(self._directory, self._basename) = os.path.split(self._directory)
			
		self._description_file = os.path.join(self._directory, self._basename, self._basename + '.desc')
		self._current_frame = start_frame
		self._end_frame = end_frame
		self._start_frame = start_frame
		self._make_current_path()	
		if not (os.path.isfile(self._coord_file) and os.path.isfile(self._voro_file)):
			raise Exception("Invalid repository: start frame out of range")
	
		self._load_balls()
		self._load_frame()
				
		if (end_frame is None):
			#if not end defined then find last accessible frame...
			#first remember start frame
			old_cur_frame = self._current_frame						
			#then try possible files
			while (os.path.isfile(self._coord_file) and os.path.isfile(self._voro_file)):
				self._current_frame += 1
				self._make_current_path()			
			#now remember it and reasign old current frame
			self._end_frame = self._current_frame - 1
			self._current_frame = old_cur_frame			
			
	def advance(self):	
		if (self._current_frame >= self._end_frame-1):	
			#we are not allowed to go further (see __init__)
			return False
		self._current_frame += 1
		self._make_current_path()
		if not (os.path.isfile(self._coord_file) and os.path.isfile(self._voro_file)):
			#undo advance...
			self._current_frame -= 1
			self._make_current_path()
			#say that we cannot go further
			return False
		
		#we have new frame so we load it
		self._load_frame()
		return True
		
	def go(self, frame):
		if frame > self._end_frame or frame < self._start_frame:
			return False
		self._current_frame = frame
		self._make_current_path()
		self._load_frame()
		return True
		
	def leap(self, offset):
		frame = self._current_frame + offset
		return self.go(frame)
		
		
	def first_frame(self):
		return self._start_frame;
		
	def last_frame(self):
		return self._end_frame;
		
	def current_frame(self):
		return self._current_frame;
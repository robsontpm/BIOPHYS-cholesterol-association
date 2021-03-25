import dmga2py
import math
from random import random
from math import pi
from math import tan
from math import sqrt

#TODO: dodac klasy reprezentujace kinetic_structure

class Geometry:
	'''
	Base (abstract) class for all geometry objects.
	Any geometry object must have bounding_box property
	that defines rectangular box that the geometry lies inside.
	bounding_box is a set of 6 float numbers 
	(min_x, min_y, min_z, max_x, max_y, max_z).
	'''
	def __init__(self):		
		self.bounding_box = (-1.0, -1.0, -1.0, 1.0, 1.0, 1.0)

	def dist(self, point_1, point_2):
		'''
		computes distance for two items in this geometry
		:param point_1: tuple/iterable (x, y, z)
		:param point_2: tuple/iterable (x, y, z)
		:return: distance between poin_1 and point_2 in this geometry
		'''
		return sqrt(sum([(point_1[k]-point_2[k])**2 for k in range(0,3)]))
		
	def transform(self, id, x, y, z, r):
		'''
		Transformation of the input coordinates.
		Usually it is identity, but some geometries (eg. cast geometries)
		may override this to do cast onto some surface.
		
		**NOTE** if writing your own Geometry you have an option to change id dynamically
		'''
		return (id, x, y, z, r)
	
	def on_boundary(self, neighbour):
		'''
		Tests if a side with a given neighbour is on boundary of the geometry
		(was created with some artificial wall not related to other particle)
		usually it is neighbour < 0.
		
		:param neighbour: a cell number, the neighbour responsible for generating some side of a cell.
		'''
		return False
	
	def adjust_area(self, area):
		'''
		Default does nothing.
		It is useful for approximate geometries, where real area may be recovered in some way
		from approximate one. 
		
		:param area: approximate area obtained from computations
		:return: adjusted area, should be exact value of the area
		'''
		return area
	
	def adjust_volume(self, volume):
		'''
		Default does nothing.
		It is useful for approximate geometries, where real volume may be recovered in some way
		from approximate one.
		
		:param volume: approximate volume obtained from computations
		:return: adjusted volume, should be exact value of the volume  
		'''
		return volume
	
class BoxGeometry(Geometry):
	'''
	This is orthogonal geometry given by bounding box
	(x_min, y_min, z_min, x_max, y_max, z_max). It may be 
	periodic in any direction.
	
	:param x_min,y_min,z_min,x_max,y_max,z_max: Float numbers
	:param x_periodic,y_periodic,z_periodic: boolean (default: False - not periodic)
	'''
	def __init__(self, x_min, y_min, z_min, x_max, y_max, z_max, x_periodic = False, y_periodic = False, z_periodic = False):
		'''
		create new instance by giving box coordinates and periodic information
		
		:param x_min,y_min,z_min,x_max,y_max,z_max: Float numbers
		:param x_periodic,y_periodic,z_periodic: boolean (default: False - not periodic)
		'''
		self.x_periodic = x_periodic 
		self.y_periodic = y_periodic
		self.z_periodic = z_periodic
		self.bounding_box = (x_min, y_min, z_min, x_max, y_max, z_max)	
		self._geometry_handle = dmga2py.new_orthogonal_geometry(x_min, y_min, z_min, x_max, y_max, z_max, x_periodic, y_periodic, z_periodic)	

	def _dist_helper(self, dist, dim_len, is_periodic):
		dist = abs(dist)
		if is_periodic and dist > 0.5 * dim_len:
			return dim_len - dist
		else:
			return dist

	def dist(self, point_1, point_2):
		'''
		computes distance for two items in this geometry. It takes into consideration
		Periodic Boundary Conditions in this Box
		:param point_1: tuple/iterable (x, y, z)
		:param point_2: tuple/iterable (x, y, z)
		:return: distance between poin_1 and point_2 in this geometry
		'''
		vals = [
			self._dist_helper(point_1[0] - point_2[0], self.bounding_box[3] - self.bounding_box[0], self.x_periodic),
			self._dist_helper(point_1[1] - point_2[1], self.bounding_box[4] - self.bounding_box[1], self.y_periodic),
			self._dist_helper(point_1[2] - point_2[2], self.bounding_box[5] - self.bounding_box[2], self.z_periodic),
		]
		return sqrt(sum(x**2 for x in vals))		
	
	def __del__(self):
		'''
		free all C++ objects if necessary
		'''
		# print "freeing geometry_handle", self._geometry_handle
		dmga2py.free_object(self._geometry_handle)	
		
	def on_boundary(self, neighbour):
		'''
		:retrun: True if neighbour id < 0, False otherwise
		''' 
		return neighbour < 0

class OrthogonalGeometry(BoxGeometry):
	'''
	Basic geometry for Molecular Dynamic simulations.
	It defines rectangular box as in PDB files (that is with x_size, y_size, z_size).
	
	The bounding box for this geometry is (0.0, 0.0, 0.0, x_size, y_size, z_size)
	
	:param x_size,y_size,z_size: Floats, size of the box
	:param x_periodic,y_periodic,z_periodic: boolean (default: False - not periodic)	
	'''
	def __init__(self, x_size, y_size, z_size, x_periodic = False, y_periodic = False, z_periodic = False):
		'''
		creates new instance
		:param x_size,y_size,z_size: Floats, size of the box
		:param x_periodic,y_periodic,z_periodic: boolean (default: False - not periodic)
		'''
		BoxGeometry.__init__(self, 0, 0, 0, x_size, y_size, z_size, x_periodic, y_periodic, z_periodic)
		self.x_size = x_size 
		self.y_size = y_size 
		self.z_size = z_size 		
		
class CastSphereGeometry(Geometry):
	'''
	This is rather articicial but still may be useful.
	It defines geometry that cast all particles on the sphere with given radius > 0
	and in given origin (default: (0,0,0)). Bounding box is computed as the minimal
	box containing given sphere. It redefines transform() to cast particles on sphere.
	it cannot be periodic (you do not have a way to define it properly though).
	
	:param radius: is double > 0
	:param origin: is tuple (default: (0,0,0))	
	'''
	def __init__(self, radius, origin = (0,0,0)):
		'''
		creates new instance
		:param radius: is double > 0
		:param origin: is tuple (default: (0,0,0))
		'''
		self.origin = origin
		self.radius = radius	
		self.bounding_box = (origin[0] - self.radius, 
							 origin[1] - self.radius, 
							 origin[2] - self.radius, 
							 origin[0] + self.radius, 
							 origin[1] + self.radius, 
							 origin[2] + self.radius)		
		self._geometry_handle = dmga2py.new_orthogonal_geometry(origin[0] - self.radius, 
																origin[1] - self.radius, 
																origin[2] - self.radius, 
																origin[0] + self.radius, 
																origin[1] + self.radius,
																origin[2] + self.radius, 
																0, 0, 0)
		self._preset_handle = dmga2py.new_geometry_preset_sphere(origin[0], origin[1], origin[2], self.radius)
		dmga2py.add_geometry_preset(self._geometry_handle, self._preset_handle)
		
	def transform(self, id, x, y, z, r):
		'''
		Cast particle to the sphere taking into account particles that lies near origin
		(those particles are not cast - so be aware)
		'''
		dist = math.sqrt( (x - self.origin[0])**2 + (y - self.origin[1])**2 + (z - self.origin[2])**2 )
 		if (dist < 0.001):
 			return (id, x, y, z, r)
# 		dist = dist + 0.02
		# scale = (1.0-random()*0.05);
		scale = 1.0
		x = scale * self.radius * (x / dist)
		y = scale * self.radius * (y / dist)
		z = scale * self.radius * (z / dist)
		return (id, x, y, z, r)	
	
	def __del__(self):
		'''
		free all C++ objects if necessary
		'''		
		# print "freeing geometry_preset_handle", self._preset_handle
		dmga2py.free_object(self._preset_handle)
		
	def on_boundary(self, neighbour):
		'''
		:retrun: True if neighbour id < 0, False otherwise
		''' 
		return neighbour < 0
		
class CastCylinderGeometry(Geometry):
	'''
	This geometry cast particles onto side of a cylinder. It may be useful in investigating
	properties of particular packaging of particles (haxagonal phase in DDPC or something similar)
	the cylinder is oriented "upwards" along z direction, and you must specify radius of the base > 0
	and its height > 0. It is by default periodic in z direction. *NOTE* that if z_periodic is set to 
	False, then particles with z > height or z < 0 will not be stored in the container!
	
	:param radius: double > 0
	:param height: double > 0
	:param z_periodic: boolean (default: True - periodic in z direction)
	:param precision: integer and defines the precision of the cylinder - higher the better but slower (default: 32) 
	'''
	def __init__(self, radius, height, z_periodic = True, precision = 32):
		'''
		creates a new instance
		:param radius: is double > 0
		:param height: is double > 0
		:param z_periodic: is boolean (default: True - periodic in z direction)
		:param precision: is integer and defines the precision of the cylinder - higher the better but slower (default: 32)
		'''
		#self.vector = vector
		#self.point = point 
		self.point = (0,0,0);
		self.vector = (0,0,1);
		self.radius = radius
		self.height = height
		self.precision = precision
# 		if (height_periodic): 
# 			height_periodic = 1
# 		else:
# 			height_periodic = 0
		#TODO: policzyc bounding box!
		#todo: zaimplementowac obrot do vector
		(x1, y1, z1, x2, y2, z2) = (-self.radius, -self.radius, 0.0, self.radius, self.radius, self.height)
		(x1, y1, z1, x2, y2, z2) = (2*x1, 2*y1, z1, 2*x2, 2*y2, z2)
		self.bounding_box = (x1, y1, z1, x2, y2, z2)
		self._geometry_handle = dmga2py.new_orthogonal_geometry(x1, y1, z1,
																x2, y2, z2,																
																False, False, z_periodic)		
#  		self._preset_handle = dmga2py.new_geometry_preset_one_cut_cylinder(
#  													self.radius,
#  													point[0], point[1], point[2], 
#  													vector[0], vector[1], vector[2])
  		self._preset_handle = dmga2py.new_geometry_preset_approximate_cylinder(self.radius, self.height, precision)
# 		self._preset_handle = dmga2py.new_geometry_preset_sphere(point[0], point[1], point[2], self.radius)
		dmga2py.add_geometry_preset(self._geometry_handle, self._preset_handle)		
		
	def transform(self, id, x, y, z, r):			
		'''
		Casts particle onto the side of the cilinder
		
		Particles close to origin (x = 0, y = 0) will not be transformed
		'''
		#do something to cast onto Cone
		origin = self.point;
		dist = math.sqrt( (x - origin[0])**2 + (y - origin[1])**2 )
 		if (dist < 0.001):
 			return (id, x, y, z, r)
		#scale = (1.0-random()*0.05);
		scale = 1.0
		x = scale * self.radius * (x / dist)
		y = scale * self.radius * (y / dist)
		if (z < 0):
			z = -z;
		return (id, x, y, z, r)	
	
	def __del__(self):
		'''
		free C++ counterpart
		''' 
		# print "freeing geometry_preset_handle", self._preset_handle
		dmga2py.free_object(self._preset_handle)
		
	def on_boundary(self, neighbour):
		'''
		:retrun: True if neighbour id < 0, False otherwise
		''' 
		return neighbour < -10 # neighbour < 0 and (neighbour < -6 or neighbour > -5)	
	
	def adjust_area(self, area):
		'''
		It can rigorously adjust computed area. 
		'''
		return area * pi / (float(self.precision) * tan(pi / float(self.precision)))
	
class CastOrthoPlaneGeometry(OrthogonalGeometry):
	'''
	This class defines a cast onto the given plane (xy), (xz) or (xz).
	It creates cast geometry around existing OrthogonalGeometry
	(It is something like Decorator)
	
	:param from_geometry: is OrthogonalGeometry
	:param vec: is one of (1,0,0), (0,1,0), (0,0,1)
	'''
	def __init__(self, from_geometry, vec):
		'''
		:param from_geometry: is OrthogonalGeometry
		:param vec: is one of (1,0,0), (0,1,0), (0,0,1)
		'''
		if (vec == (1,0,0)):
			invvec = (0,1,1)
			self._on_boundary = -1
		elif (vec == (0,1,0)):
			invvec = (1,0,1)
			self._on_boundary = -3
		elif (vec == (0,0,1)):
			invvec = (1,1,0)
			self._on_boundary = -5
		else:
			raise Exception("CastOrthoPlaneGeometry.__init__(): Bad cast vector!")
		x_size = from_geometry.x_size * invvec[0] + vec[0];
		y_size = from_geometry.y_size * invvec[1] + vec[1];
		z_size = from_geometry.z_size * invvec[2] + vec[2];
		x_periodic = from_geometry.x_periodic and invvec[0]
		y_periodic = from_geometry.y_periodic and invvec[1]
		z_periodic = from_geometry.z_periodic and invvec[2] 
		OrthogonalGeometry.__init__(self, x_size, y_size, z_size, x_periodic, y_periodic, z_periodic)
		self.vector = vec
		self.cast_vector = invvec;		
		
	def transform(self, id, x, y, z, r):
		'''
		cast particle to a defined plane
		'''
		return (id, x * self.cast_vector[0], y * self.cast_vector[1], z * self.cast_vector[2], r)
	
	def on_boundary(self, neighbour):
		'''
		may be used to filter only geometry on the casting plane
		'''
		return neighbour == self._on_boundary
		
class CastPlaneGeometry(OrthogonalGeometry):
	'''
	NotImplementedYet due to some problems with definition of the cast
	'''
	def __init__(self, from_geometry, cast_vector, x_periodic = False, y_periodic = False):
		raise Exception("CastPlaneGeometry.__init__(): NotImplementedYet")
		self.from_geomtery = from_geometry
		self.cast_vector = cast_vector
		self.x_periodic = x_periodic
		self.y_periodic = y_periodic
		#TODO: policzyc bounding box!
		self.bounding_box = (0, 0, 0, 0, 0, 0)
		
	def transform(self, id, x, y, z, r):
		#do something to cast onto Plane
		return (id, x, y, z, r)					

	

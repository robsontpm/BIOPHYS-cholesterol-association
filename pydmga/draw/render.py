import random
from math import pi, sin, cos, sqrt, asin
import math
import os #os.path.*
import re #regular expressions (split)
import sys #sys.argv

from panda3d.core import * 
from pandac.PandaModules import NodePath
from pandac.PandaModules import Point3
from pandac.PandaModules import Vec4
from pandac.PandaModules import LineSegs
from direct.showbase.ShowBase import ShowBase
from direct.task import Task

from geometry import * #voro.geometry
import assets #voro.graphics

class Color:
	'''
	predefined colors to use with render
	'''
	GREEN_OPAQUE   = (0.0,1.0,0.0,1.0)
	GREEN10_OPAQUE = (0.0,0.1,0.0,1.0)
	GREEN20_OPAQUE = (0.0,0.2,0.0,1.0)
	GREEN30_OPAQUE = (0.0,0.3,0.0,1.0)
	GREEN40_OPAQUE = (0.0,0.4,0.0,1.0)
	GREEN50_OPAQUE = (0.0,0.5,0.0,1.0)
	GREEN60_OPAQUE = (0.0,0.6,0.0,1.0)
	GREEN70_OPAQUE = (0.0,0.7,0.0,1.0)
	GREEN80_OPAQUE = (0.0,0.8,0.0,1.0)
	GREEN90_OPAQUE = (0.0,0.9,0.0,1.0)
	
	RED_OPAQUE   = (1.0,0.0,0.0,1.0)
	RED10_OPAQUE = (0.1,0.0,0.0,1.0)
	RED20_OPAQUE = (0.2,0.0,0.0,1.0)
	RED30_OPAQUE = (0.3,0.0,0.0,1.0)
	RED40_OPAQUE = (0.4,0.0,0.0,1.0)
	RED50_OPAQUE = (0.5,0.0,0.0,1.0)
	RED60_OPAQUE = (0.6,0.0,0.0,1.0)
	RED70_OPAQUE = (0.7,0.0,0.0,1.0)
	RED80_OPAQUE = (0.8,0.0,0.0,1.0)
	RED90_OPAQUE = (0.9,0.0,0.0,1.0)
	
	BLUE_OPAQUE   = (0.0,0.0,1.0,1.0)
	BLUE10_OPAQUE = (0.0,0.0,0.1,1.0)
	BLUE20_OPAQUE = (0.0,0.0,0.2,1.0)
	BLUE30_OPAQUE = (0.0,0.0,0.3,1.0)
	BLUE40_OPAQUE = (0.0,0.0,0.4,1.0)
	BLUE50_OPAQUE = (0.0,0.0,0.5,1.0)
	BLUE60_OPAQUE = (0.0,0.0,0.6,1.0)
	BLUE70_OPAQUE = (0.0,0.0,0.7,1.0)
	BLUE80_OPAQUE = (0.0,0.0,0.8,1.0)
	BLUE90_OPAQUE = (0.0,0.0,0.9,1.0)
	
	WHITE_OPAQUE  = (1.0,1.0,1.0,1.0)
	GRAY10_OPAQUE = (0.1,0.1,0.1,1.0)
	GRAY20_OPAQUE = (0.2,0.2,0.2,1.0)
	GRAY30_OPAQUE = (0.3,0.3,0.3,1.0)
	GRAY40_OPAQUE = (0.4,0.4,0.4,1.0)
	GRAY50_OPAQUE = (0.5,0.5,0.5,1.0)
	GRAY60_OPAQUE = (0.6,0.6,0.6,1.0)
	GRAY70_OPAQUE = (0.7,0.7,0.7,1.0)
	GRAY80_OPAQUE = (0.8,0.8,0.8,1.0)
	GRAY90_OPAQUE = (0.9,0.9,0.9,1.0)	
	BLACK_OPAQUE  = (0.0,0.0,0.0,1.0)
	
	GREEN_HALF   = (0.0,1.0,0.0,0.5)
	GREEN10_HALF = (0.0,0.1,0.0,0.5)
	GREEN20_HALF = (0.0,0.2,0.0,0.5)
	GREEN30_HALF = (0.0,0.3,0.0,0.5)
	GREEN40_HALF = (0.0,0.4,0.0,0.5)
	GREEN50_HALF = (0.0,0.5,0.0,0.5)
	GREEN60_HALF = (0.0,0.6,0.0,0.5)
	GREEN70_HALF = (0.0,0.7,0.0,0.5)
	GREEN80_HALF = (0.0,0.8,0.0,0.5)
	GREEN90_HALF = (0.0,0.9,0.0,0.5)
	
	RED_HALF   = (1.0,0.0,0.0,0.5)
	RED10_HALF = (0.1,0.0,0.0,0.5)
	RED20_HALF = (0.2,0.0,0.0,0.5)
	RED30_HALF = (0.3,0.0,0.0,0.5)
	RED40_HALF = (0.4,0.0,0.0,0.5)
	RED50_HALF = (0.5,0.0,0.0,0.5)
	RED60_HALF = (0.6,0.0,0.0,0.5)
	RED70_HALF = (0.7,0.0,0.0,0.5)
	RED80_HALF = (0.8,0.0,0.0,0.5)
	RED90_HALF = (0.9,0.0,0.0,0.5)
	
	BLUE_HALF   = (0.0,0.0,1.0,0.5)
	BLUE10_HALF = (0.0,0.0,0.1,0.5)
	BLUE20_HALF = (0.0,0.0,0.2,0.5)
	BLUE30_HALF = (0.0,0.0,0.3,0.5)
	BLUE40_HALF = (0.0,0.0,0.4,0.5)
	BLUE50_HALF = (0.0,0.0,0.5,0.5)
	BLUE60_HALF = (0.0,0.0,0.6,0.5)
	BLUE70_HALF = (0.0,0.0,0.7,0.5)
	BLUE80_HALF = (0.0,0.0,0.8,0.5)
	BLUE90_HALF = (0.0,0.0,0.9,0.5)
	
	WHITE_HALF  = (1.0,1.0,1.0,0.5)
	GRAY10_HALF = (0.1,0.1,0.1,0.5)
	GRAY20_HALF = (0.2,0.2,0.2,0.5)
	GRAY30_HALF = (0.3,0.3,0.3,0.5)
	GRAY40_HALF = (0.4,0.4,0.4,0.5)
	GRAY50_HALF = (0.5,0.5,0.5,0.5)
	GRAY60_HALF = (0.6,0.6,0.6,0.5)
	GRAY70_HALF = (0.7,0.7,0.7,0.5)
	GRAY80_HALF = (0.8,0.8,0.8,0.5)
	GRAY90_HALF = (0.9,0.9,0.9,0.5)	
	BLACK_HALF  = (0.0,0.0,0.0,0.5)
	
	GREEN_QUARTER   = (0.0,1.0,0.0,0.25)
	GREEN10_QUARTER = (0.0,0.1,0.0,0.25)
	GREEN20_QUARTER = (0.0,0.2,0.0,0.25)
	GREEN30_QUARTER = (0.0,0.3,0.0,0.25)
	GREEN40_QUARTER = (0.0,0.4,0.0,0.25)
	GREEN50_QUARTER = (0.0,0.5,0.0,0.25)
	GREEN60_QUARTER = (0.0,0.6,0.0,0.25)
	GREEN70_QUARTER = (0.0,0.7,0.0,0.25)
	GREEN80_QUARTER = (0.0,0.8,0.0,0.25)
	GREEN90_QUARTER = (0.0,0.9,0.0,0.25)
	
	RED_QUARTER   = (1.0,0.0,0.0,0.25)
	RED10_QUARTER = (0.1,0.0,0.0,0.25)
	RED20_QUARTER = (0.2,0.0,0.0,0.25)
	RED30_QUARTER = (0.3,0.0,0.0,0.25)
	RED40_QUARTER = (0.4,0.0,0.0,0.25)
	RED50_QUARTER = (0.5,0.0,0.0,0.25)
	RED60_QUARTER = (0.6,0.0,0.0,0.25)
	RED70_QUARTER = (0.7,0.0,0.0,0.25)
	RED80_QUARTER = (0.8,0.0,0.0,0.25)
	RED90_QUARTER = (0.9,0.0,0.0,0.25)
	
	BLUE_QUARTER   = (0.0,0.0,1.0,0.25)
	BLUE10_QUARTER = (0.0,0.0,0.1,0.25)
	BLUE20_QUARTER = (0.0,0.0,0.2,0.25)
	BLUE30_QUARTER = (0.0,0.0,0.3,0.25)
	BLUE40_QUARTER = (0.0,0.0,0.4,0.25)
	BLUE50_QUARTER = (0.0,0.0,0.5,0.25)
	BLUE60_QUARTER = (0.0,0.0,0.6,0.25)
	BLUE70_QUARTER = (0.0,0.0,0.7,0.25)
	BLUE80_QUARTER = (0.0,0.0,0.8,0.25)
	BLUE90_QUARTER = (0.0,0.0,0.9,0.25)
	
	WHITE_QUARTER  = (1.0,1.0,1.0,0.25)
	GRAY10_QUARTER = (0.1,0.1,0.1,0.25)
	GRAY20_QUARTER = (0.2,0.2,0.2,0.25)
	GRAY30_QUARTER = (0.3,0.3,0.3,0.25)
	GRAY40_QUARTER = (0.4,0.4,0.4,0.25)
	GRAY50_QUARTER = (0.5,0.5,0.5,0.25)
	GRAY60_QUARTER = (0.6,0.6,0.6,0.25)
	GRAY70_QUARTER = (0.7,0.7,0.7,0.25)
	GRAY80_QUARTER = (0.8,0.8,0.8,0.25)
	GRAY90_QUARTER = (0.9,0.9,0.9,0.25)	
	BLACK_QUARTER  = (0.0,0.0,0.0,0.25)	
	

class Primitive(object):
	'''basic class for all primitives that Render class can handle'''	
	def attach_to(self, render):
		'''each geometry object should overwrite this'''
		pass

class PointCloud(Primitive):
	'''This can be used to render a set of points'''
	
	def __init__(self):
		'''creates empty set'''
		format_array = GeomVertexArrayFormat()
		format_array.addColumn(InternalName.make('vertex'), 3, Geom.NTFloat32, Geom.CPoint)
		format_array.addColumn(InternalName.make('size'), 1, Geom.NTFloat32, Geom.COther)
		format_array.addColumn(InternalName.make('color'), 4, Geom.NTFloat32, Geom.CColor)	
		format = GeomVertexFormat()
		format.addArray(format_array)				
		format = GeomVertexFormat.registerFormat(format)
		self.vertices_raw = GeomVertexData('vertices', format, Geom.UHStatic);			
		self.vertex_writer = GeomVertexWriter(self.vertices_raw, 'vertex')
		self.size_writer = GeomVertexWriter(self.vertices_raw, 'size')
		self.color_writer = GeomVertexWriter(self.vertices_raw, 'color')
		self.vertices_primitive = GeomPoints(Geom.UHStatic)	
		self.vertex_count = 0
		
	def add_vertex(self, point, size = 1.0, color = None):		
		'''adds new vertex with given color and size s'''
		if color is None: color = (1.0, 1.0, 1.0, 1.0)
		(r,g,b,a) = color;
		(x,y,z) = point
		self.vertices_primitive.addVertex(self.vertex_count)
		self.vertex_count += 1
		self.vertex_writer.addData3f(x, y, z)
		self.size_writer.addData1f(size)
		self.color_writer.addData4f(r, g, b, a)
		return self
		
	def add_vertices(self, points, size, color):
		'''adds many points all in the same color'''
		for p in points:
			self.add_vertex(p, size, color)
		return self
		
	def attach_to(self, render):
		'''creates Panda3D representation that can be attached to render
		this function is mainly used inside Render class and its derivatives and
		user not need to worry to call it'''
		self.vertices_primitive.closePrimitive()				
		vertices_geom = Geom(self.vertices_raw)
		vertices_geom.addPrimitive(self.vertices_primitive) 
		vertices_node = GeomNode('vertices_node')
		vertices_node.addGeom(vertices_geom) 	
		node = render.render.attachNewNode(vertices_node)
		return node
		
class ConvexPolygons(Primitive):
	'''It allows to create a set of polygons with common points,
	for example it works perfectly for Voronoi Diagram Cells'''
	
	def __init__(self, points = None, color = None):
		'''
		you can skip points and color
		if you pass array of points then the constructor automatically
		creates one polygon on all points in a given order
		if color is skipped then white non-translucent color is used
		'''
		format_array = GeomVertexArrayFormat()
		format_array.addColumn(InternalName.make('vertex'), 3, Geom.NTFloat32, Geom.CPoint)
		format_array.addColumn(InternalName.make('color'), 4, Geom.NTFloat32, Geom.CColor)	
		format = GeomVertexFormat()
		format.addArray(format_array)				
		format = GeomVertexFormat.registerFormat(format)
		self.vertices_raw = GeomVertexData('polygon', format, Geom.UHStatic);			
		self.vertex_writer = GeomVertexWriter(self.vertices_raw, 'vertex')
		self.color_writer = GeomVertexWriter(self.vertices_raw, 'color')
		self.triangles_primitive = GeomTriangles(Geom.UHStatic)	
		self.vertex_count = 0
		if (color is None):
			color = (1.0, 1.0, 1.0, 1.0)
			self.base_color = color
		self.start_side()
		if (points is not None):
			for i, p in enumerate(points):
				self.add_vertex(p, color)
				self.add_side_point(i)
		
	def add_vertex(self, point, color = None):
		'''stores new vertex with given color'''
		if color is None: color = self.base_color
		(r,g,b,a) = color
		(x,y,z) = point
		self.vertex_count += 1
		self.vertex_writer.addData3f(x, y, z)
		self.color_writer.addData4f(r, g, b, a)
		return self
		
	def add_vertices(self, vertices, color = None):
		'''adds many points all in the same color'''
		for p in vertices:
			self.add_vertex(p, color)
		return self
	
	def add_side_points(self, numbers):
		'''adds all vertices with numbers in number to a current side'''
		for n in numbers:
			self.add_side_point(number)
		return self
		
	def add_side(self, vertices, color = None):
		'''creates a whole new side with all new vertices'''
		curr = self.vertex_count
		self.start_side()
		for p in vertices:
			self.add_vertex(p, color)
			self.add_side_point(curr)
			curr += 1
		return self		
		
	def start_side(self):
		'''closes last side and starts new one (new polygon), 
		all subsequent calls to add_side_point will be addressing this 
		newly created side'''
		self.base_point = None
		self.second_point = None
		return self		
		
	def add_side_point(self, number):
		'''adds new point to a side'''
		if (self.base_point is None):
			self.base_point = number
		elif (self.second_point is None):
			self.second_point = number
		else:
			self.triangles_primitive.addVertices(self.base_point, self.second_point, number)
			self.second_point = number;
		return self
		
	def attach_to(self, render):
		'''creates Panda3D representation that can be attached to render
		this function is mainly used inside Render class and its derivatives and
		user not need to worry to call it'''
		self.triangles_primitive.closePrimitive()				
		poly_geom = Geom(self.vertices_raw)
		poly_geom.addPrimitive(self.triangles_primitive) 
		poly_node = GeomNode('poly_node')
		poly_node.addGeom(poly_geom) 
		node = render.render.attachNewNode(poly_node)
		node.setTransparency(TransparencyAttrib.MDual, 1)
		return node
	
class Outline(Primitive):
	def __init__(self, points, thickness = 1.0, color = None):
		'''creates a closed contour with given points'''
		self.points = points
		self.thickness = thickness
		if (color is None): color = (0.5, 0.5, 0.5, 1.0)
		self.color = color
		
	def attach_to(self, render):
		'''creates Panda3D representation that can be attached to render
		this function is mainly used inside Render class and its derivatives and
		user not need to worry to call it'''
		new_line_seqs = LineSegs()
		new_line_seqs.setThickness(self.thickness)		
		(r,g,b,a) = self.color
		new_line_seqs.setColor(Vec4(r,g,b,a))			
		(x,y,z) = self.points[0]
		new_line_seqs.moveTo(Point3(x, y, z))
		for i in xrange(1,len(self.points)):
			(x,y,z) = self.points[i]
			new_line_seqs.drawTo(Point3(x, y, z))
		(x,y,z) = self.points[0]
		new_line_seqs.drawTo(Point3(x, y, z))
		lines = NodePath(new_line_seqs.create())
		lines.reparentTo(render.render)
		return lines
	
class Lines(Primitive):
	def __init__(self, points, thickness = 1.0, color = None):
		'''creates a poly line with given points'''
		self.points = points
		self.thickness = thickness
		if (color is None): color = (0.5, 0.5, 0.5, 1.0)
		self.color = color
		
	def attach_to(self, render):
		'''creates Panda3D representation that can be attached to render
		this function is mainly used inside Render class and its derivatives and
		user not need to worry to call it'''
		new_line_seqs = LineSegs()
		new_line_seqs.setThickness(self.thickness)		
		(r,g,b,a) = self.color
		new_line_seqs.setColor(Vec4(r,g,b,a))			
		(x,y,z) = self.points[0]
		new_line_seqs.moveTo(Point3(x, y, z))
		for i in xrange(1,len(self.points)):
			(x,y,z) = self.points[i]
			new_line_seqs.drawTo(Point3(x, y, z))
		lines = NodePath(new_line_seqs.create())
		lines.reparentTo(render.render)
		return lines
	
		
class Asset(Primitive):
	'''this represents a class of all object loaded from external files'''
	def __init__(self, loader, point, scale, color = None):
		self.x = point[0]
		self.y = point[1]
		self.z = point[2]
		self.r = scale
		self.color = color
		self.loader = loader
		
	def attach_to(self, render):
		'''creates Panda3D representation that can be attached to render
		this function is mainly used inside Render class and its derivatives and
		user not need to worry to call it'''
		model = self.loader(render.loader, self.x, self.y, self.z, self.r)
		if not(self.color is None):
			(r,g,b,a) = self.color 
			model.setColor(r,g,b,a)  			
		model.reparentTo(render.render)
		model.setTransparency(TransparencyAttrib.MDual, 1)
		return model
	
class OnSpherePrimitive(Primitive):
	def __init__(self, origin, radius):
		self.origin = origin
		self.radius = radius
		
	def cast_to_sphere(self, vertex):
		v = [(vertex[i], self.origin[i]) for i in xrange(3)] # make a list of pairs suitable for computations with fancy python functions ;)
		d = sqrt(sum(a for a in map(lambda(x):(x[0]-x[1])**2, v))) # compute length of vector
		result = map(lambda(x): (x[0]-x[1]) / d * self.radius + x[1], v)
		return result # cast to sphere	
	
	def distance(self, v, u = (0,0,0)):
		return sqrt(sum((v[i] - u[i])**2 for i in range(3)))
	
	def normalize(self, v, u = (0,0,0)):
		d = self.distance(v,u)
		return [(v[i] - u[i]) / d for i in range(3)]
	
	def dot(self, u, v, c = (0,0,0)):
		return sum((c[i] - u[i]) * (v[i] - c[i]) for i in range(3))
		
	def cross(self, u, v, c = (0,0,0)):
		return (
		(u[1]-c[1])*(v[2]-c[2]) - (u[2]-c[2])*(v[1]-c[1]),
		(u[2]-c[2])*(v[0]-c[0]) - (u[0]-c[0])*(v[2]-c[2]),
		(u[0]-c[0])*(v[1]-c[1]) - (u[1]-c[1])*(v[0]-c[0]))
		
# 		(c[1]-u[1])*(v[2]-c[2]) - (c[2]-u[2])*(v[1]-c[1]),
# 		(c[2]-u[2])*(v[0]-c[0]) - (c[0]-u[0])*(v[2]-c[2]),
# 		(c[0]-u[0])*(v[1]-c[1]) - (c[1]-u[1])*(v[1]-c[1]))		
		
	
class SphereArcs(OnSpherePrimitive):
	'''this allow to draw part of circles that are the result of intersection of two spheres
	to define this you need to give three points that lie on the plane defined by two intersecting 
	balls 
	'''
	def __init__(self, origin, radius, thickness = 1.0, color = (0.0, 1.0, 0.0, 1.0)):
		super(SphereArcs, self).__init__(origin, radius)		
		self.thickness = thickness
		self.color = color
		self.vertex_count = 0
		self.vertices = []
		self.line_seqs = LineSegs()
		self.line_seqs.setThickness(self.thickness)
		(r,g,b,a) = self.color
		self.line_seqs.setColor(Vec4(r,g,b,a))
		
	def _generate_points_on_arc(self, center, from_v, to_v, circle_r = None):
		'''we assume that the arc is centered at center, and is going from
		from_v to to_v in a given direction (for example counterclockwise - to set)''' 
		dist_sq = sum( (from_v[i] - to_v[i])**2 for i in xrange(3)) # distance between points
		if circle_r is None:
			circle_r = sqrt(sum((from_v[i] - center[i])**2 for i in range(3)))		
		if (dist_sq < (self.radius / 10.0)**2):
			return [to_v]
		else:
			new_v = [ ((from_v[i] + to_v[i]) / 2.0 - center[i]) for i in xrange(3)]
			d = sqrt(sum(new_v[i]**2 for i in range(3)))
			new_v = [center[i] + (new_v[i] / d * circle_r) for i in range(3)]			
			a = self._generate_points_on_arc(center, from_v, new_v, circle_r)
			b = self._generate_points_on_arc(center, new_v, to_v, circle_r)
			return a + b
	
	def _is_greater_arc(self, center, from_v, to_v):
		return False #TODO: stub!!!!!!
	
	def _make_arc(self, center, from_v, to_v):
		'''
		if self._is_greater_arc(center, from_v, to_v):
			v_prim = self.cast_to_sphere([3.0 * self.origin[i] - from_v[i] + to_v[i] for i in range(3)]) #some point on the greater arc between from and to (it is computed using: cat_to_sphere(origin + 2 * (origin - mid_point(from, to))   (origin - mid_point(from, to)) is a vector from midpoint to origin 
			list = self._generate_points_on_arc(center, from_v, v_prim)
			list += self._generate_points_on_arc(center, v_prim, to_v)
		else:
			list = self._generate_points_on_arc(center, from_v, to_v)
		(x,y,z) = from_v
		self.line_seqs.moveTo(Point3(x, y, z))
		for p in list:
			(x,y,z) = p
			self.line_seqs.drawTo(Point3(x, y, z))
		'''
		r = self.distance(center, from_v)
		a = self.normalize(from_v, center)
		b = self.normalize(to_v, center)
		divisions = max(int(self.radius * 10), 6)
		base_theta = (pi / (2.0 * divisions))
		c = center
		(x,y,z) = from_v
		self.line_seqs.moveTo(Point3(x, y, z))
		for i in xrange(1, divisions+1):
			theta = base_theta * float(i)
			(x,y,z) = [c[i] + r * cos(theta) * a[i] + r * sin(theta) * b[i] for i in range(3)]
			self.line_seqs.drawTo(Point3(x, y, z))			
	
	def add_arc(self, center, from_v, to_v):
		'''center = of the circle, the cutting plane
		from_v = first vertex on arc, should lie on the cutting plane and on the sphere
		to_t - last vertex on arc, should lie on the cutting plane and on the sphere
		moreover we require that dist(center, from_v) == dist(center, to_v)
		'''
		#self._make_arc(center, from_v, to_v)
		#normal = (center[0]-self.origin[0], center[1]-self.origin[1], center[2]-self.origin[2])
		normal = self.normalize(center, self.origin)
		dist = self.distance(center, self.origin)
		self.add_cutting_plane_arc(normal, dist, from_v, to_v)
		return self
	
	def add_cutting_plane_arc(self, normal, distance, from_v, to_v):
		'''giving other center and distance gives us chance to handle degenerate case when center = origin of sphere
		'''
		n = self.normalize(normal) # assure normal is normal
		c = [self.origin[i] + n[i] * distance for i in range(3)] # the center of the arc
		#TODO: rzutowac from_v i to_v na okrag przeciecia plaszczyzna - ale to skomplikowane
		
		divisions = 32 # max(int(self.radius * 10), 16)
		base_theta = 2.0 * pi / divisions
		
		r = self.distance(c, from_v) 
		u = self.normalize(from_v, c)
		v = self.normalize(to_v, c)
		uxv = self.cross(u, v)		
		sin_theta = self.distance(uxv)
		cos_theta = self.dot(u, v)
		sin_sign = self.dot(n, uxv) # this should be 1, 0 or -1, 0 not matter (then sin = 0, parallel things,
		a = u
		b = self.cross(n,u)
		self.normalize(a)
		self.normalize(b)
		b = [-b[i] for i in range(3)]
		
		if (from_v == to_v):
			final_theta = 2.0 * pi
		else:
			if (cos_theta * sin_sign > 0):
				final_theta = asin(1.0 - sin_theta)
			else:
				final_theta = asin(sin_theta)
			print "add_cutting_plane_arc", final_theta, ",", sin_theta
	 		if (sin_sign < 0): #this means that we are on the other side of the line and the angle is 180 deegres plus angle complementar to asin(theta)
	 			#final_theta = asin(1.0 - sin_theta)
			 	final_theta += pi
			if (sin_sign == 0 and cos_theta > 0):
				final_theta += pi
	 		if (cos_theta * sin_sign > 0):
	 		 	final_theta += (pi / 2.0)
 		 	 		
		(x,y,z) = from_v
		self.line_seqs.moveTo(Point3(x, y, z))
		theta = base_theta
		print "add_cutting_plane_arc", "\n\tu:", u, "\n\tv:", v, "\n\tuxv:", uxv, "\n\tsin_theta:", sin_theta, "\n\tcos_theta:", cos_theta, "\n\tsign_theta:", sin_sign, "\n\ta:", a, "\n\tb:", b, "\n\tfinal_theta:", final_theta
		while theta < final_theta:
			(x,y,z) = [c[i] + r * cos(theta) * a[i] + r * sin(theta) * b[i] for i in range(3)]
			self.line_seqs.drawTo(Point3(x, y, z))
			theta += base_theta
		(x,y,z) = to_v
		self.line_seqs.drawTo(Point3(x, y, z))
		
# 		self.line_seqs.moveTo(Point3(from_v[0], from_v[1], from_v[2]))
# 		self.line_seqs.drawTo(Point3(c[0], c[1], c[2]))
# 		self.line_seqs.drawTo(Point3(to_v[0], to_v[1], to_v[2])) 
# 		if (sin_theta == 0):
# 			if (from_v != to_v):
# 				while theta <= pi:
# 					(x,y,z) = [c[i] + r * cos(theta) * a[i] + r * sin(theta) * b[i] for i in range(3)]
# 					self.line_seqs.drawTo(Point3(x, y, z))
# 					theta += base_theta
# 				(x,y,z) = to_v
# 				self.line_seqs.drawTo(Point3(x, y, z))
# 		elif (sin_theta < 0):
# 			#do half an arc			
# 			while theta <= pi:
# 				(x,y,z) = [c[i] + r * cos(theta) * a[i] + r * sin(theta) * b[i] for i in range(3)]
# 				self.line_seqs.drawTo(Point3(x, y, z))
# 				theta += base_theta
# 			#do the rest
# 			count = 0
#  			#while sin(theta) > sin_theta and count < 40: 
#  			while sin(theta) > sin_theta:
#  				(x,y,z) = [c[i] + r * cos(theta) * a[i] + r * sin(theta) * b[i] for i in range(3)]
#  				self.line_seqs.drawTo(Point3(x, y, z))
#  				theta += base_theta
#  				count += 1
#  			(x,y,z) = to_v
# 			self.line_seqs.drawTo(Point3(x, y, z))
# 		else:
# 			count = 0
# 			#while sin(theta) < sin_theta and count < 40:
# 			while sin(theta) < sin_theta:
# 				(x,y,z) = [c[i] + r * cos(theta) * a[i] + r * sin(theta) * b[i] for i in range(3)]
# 				self.line_seqs.drawTo(Point3(x, y, z))
# 				theta += base_theta
# 				count+=1
# 			(x,y,z) = to_v
# 			self.line_seqs.drawTo(Point3(x, y, z))
		return self
				
	
	def attach_to(self, render):
		'''creates Panda3D representation that can be attached to render
		this function is mainly used inside Render class and its derivatives and
		user not need to worry to call it'''
		lines = NodePath(self.line_seqs.create())
		lines.reparentTo(render.render)
		return lines	
	
	
class GreatArcs(OnSpherePrimitive):
	'''this class can draw a collection of great arcs on the given sphere
	(that is a part of the biggest possible circle on a sphere, which has origin at the origin of a sphere)'''
	def __init__(self, origin, radius, thickness = 1.0, color = (0.0, 1.0, 0.0, 1.0)):
		super(GreatArcs, self).__init__(origin, radius)
		self.thickness = thickness
		self.color = color
		self.vertex_count = 0
		self.vertices = []
		self.line_seqs = LineSegs()
		self.line_seqs.setThickness(self.thickness)
		(r,g,b,a) = self.color
		self.line_seqs.setColor(Vec4(r,g,b,a))
		
	def _generate_points_on_arc(self, from_v, to_v, depth = 0):
		'''from_v and to_v are vectors, i need to generate
		approximate middle points on surface of the sphare and add them to 
		line_seqs
		we assume that from_v, to_v are on the sphere... 
		at the end we return the list without FIRST_VERTEX (from_v), it is
		convenient in _make_arc, as we neet do do one  move, and then we can lineTo
		to all elements in this list
		'''		
		#there is a problem when from_v and to_v are antipodal... TODO: poprawic?
		dist_sq = sum( (from_v[i] - to_v[i])**2 for i in xrange(3)) # distance between points
		if (dist_sq < (self.radius / 5.0)**2 or depth > 10): # its quesss...			
			# simply connects
			return [to_v]
		else:
			# generate middle point and call recursively
			new_v = self.cast_to_sphere([ (from_v[i] + to_v[i]) / 2.0 for i in xrange(3)])
			a = self._generate_points_on_arc(from_v, new_v, depth + 1)
			b = self._generate_points_on_arc(new_v, to_v, depth + 1)
			return a + b # concat lists...  
		
	def _make_arc(self, from_v, to_v):		
		list = self._generate_points_on_arc(from_v, to_v)
		(x,y,z) = from_v
		self.line_seqs.moveTo(Point3(x, y, z))
		for p in list:
			(x,y,z) = p
			self.line_seqs.drawTo(Point3(x, y, z))

	def _make_arcs(self, vertices):
		'''TODO: przetestowac, bo narazie napisalem...'''
		first = vertices[0]
		current = first
		list = []
		for next in vertices[1:]:
			list += self._generate_points_on_arc(current, next)
			current = next		
		(x,y,z) = first
		self.line_seqs.moveTo(Point3(x, y, z))
		for p in list:
			(x,y,z) = p
			self.line_seqs.drawTo(Point3(x, y, z))			
		
	def add_vertex(self, vertex, cast_to_sphere = False):
		'''adds a single point to this object
		you can connect them wit connect() function
		'''
		self.vertex_count += 1
		if (cast_to_sphere):
			self.vertices.append(self.cast_to_sphere(vertex))
		else:
			self.vertices.append(vertex)
		return self
	
	def add_vertices(self, vertices, cast_to_sphere = False):
		'''adds a list of vertices
		you can connect them wit connect() function
		'''
		for p in vertices:
			self.add_vertex(p, cast_to_sphere)
		return self
	
	def connect(self, indices_list, closed = True):
		'''connects given list of indices by set of fragments of great arcs
		the default behaviour is that it creates a closed curve
		'''
		current = indices[0]
		for next in indices_list[1:]:
			self._make_arc(self.vertices[current], self.vertices[next])
			current = next
		if closed:
			self._make_arc(self.vertices[current], self.vertices[indices[0]])
		return self
	
	def add_arcs(self, vertices, closed = True, cast_to_sphere = False):
		'''adds a collection of fragments of great arcs between given points
		points should be on the surface of the ball
		the default behaviour is that it creates a closed curve
		'''
		first = self.vertex_count
		current = first
		self.add_vertex(vertices[0], cast_to_sphere)
		for p in vertices[1:]:
			next = self.vertex_count
			self.add_vertex(p, cast_to_sphere)
			self._make_arc(self.vertices[current], self.vertices[next])
			current = next
		if closed:
			self._make_arc(self.vertices[current], self.vertices[first])
		return self
	
	def attach_to(self, render):
		'''creates Panda3D representation that can be attached to render
		this function is mainly used inside Render class and its derivatives and
		user not need to worry to call it'''
		lines = NodePath(self.line_seqs.create())
		lines.reparentTo(render.render)
		return lines
		
class Render(ShowBase):
	'''simple class for rendering things using Panda3D'''
	def __init__(self, bgcolor = None, ltcolor = None):
		ShowBase.__init__(self)
		if bgcolor is None: bgcolor = (0.0, 0.0, 0.0, 1.0) #default is black
		if ltcolor is None: ltcolor = (1.0, 1.0, 1.0, 1.0) #default if white
		self.current_nodes = []	
		self.alight = AmbientLight('alight')		
		(r,g,b,a) = ltcolor
		self.alight.setColor(VBase4(r,g,b,a))
		self.alnp = render.attachNewNode(self.alight)
		render.setLight(self.alnp)
		(r,g,b,a) = bgcolor
		base.setBackgroundColor(r,g,b,a)
		
	def add(self, primitive):
		'''adds a class of type Primitive to this render object'''
		node = primitive.attach_to(self)
		self.current_nodes.append(node)		
		
	def clear(self):
		'''removes all nodes from render'''
		for node in self.current_nodes:
			node.removeNode()
		self.current_nodes = []


class VoroRender(Render):
 
	def __init__(self, input):
		Render.__init__(self)	
		
		self.input = input	
				
		#self.oobe()			
		self.taskMgr.add(self.frame_iterator_task, "frame_render")
		self.is_running = False
		self.need_update = True
		self.alpha = 0.0
		self.grow_particles = False
		self.show_shape = False
		self.show_proximity = False
		self.selected_cells = []
		
		self.observers = []
		
		self.balls = []	
		self.lines = None
		self.lines_thin = None
		self.lines_thick = None
		self.vertices = None
		self.sides = None
		#self._voro_draw()					
		
		for i, ball in enumerate(self.input.balls):
			model = assets.sphere(self.loader, ball.x, ball.y, ball.z, ball.r) 			
			model.reparentTo(self.render)
			self.balls.append(model)
			self.selected_cells.append(True)
		#print self.selected_cells
			
	def _voro_draw(self):
		new_line_seqs = LineSegs()
		new_line_seqs.setThickness(2.0) #TODO
		new_line_seqs.setColor(Vec4(1,1,0,0.5)) #TODO			
		n = 0
		for (points, faces) in self.input.current.voro:			
			if self.selected_cells[n]:				
				for (poly, poly_data) in faces:
					(x,y,z,point_data) = points[poly[0][0]]
					new_line_seqs.moveTo(Point3(x, y, z))
					for i in xrange(1,len(poly)):
						(x,y,z,point_data) = points[poly[i][0]]
						new_line_seqs.drawTo(Point3(x, y, z))
					(x,y,z,point_data) = points[poly[0][0]]
					new_line_seqs.drawTo(Point3(x, y, z))
				if (self.grow_particles):
					self.balls[n].setColor(1.0, 1.0, 1.0, 0.25)	
				else:
					self.balls[n].setColor(1.0, 1.0, 1.0, 1.0)	
			else:
				if (self.grow_particles):
					self.balls[n].setColor(1.0, 1.0, 1.0, 0.05)	
				else:
					self.balls[n].setColor(1.0, 1.0, 1.0, 0.5)	
			n = n + 1
		
		if not(self.vertices is None):
			self.vertices.removeNode()
			self.vertices = None		
		if not(self.sides is None):
			self.sides.removeNode()
			self.sides = None			
		if not(self.lines_thin is None):
			self.lines_thin.removeNode()
			self.lines_thin = None
		if not(self.lines_thick is None):
			self.lines_thick.removeNode()		
			self.lines_thick = None
		if not(self.lines is None):
			self.lines.removeNode()
		self.lines = NodePath(new_line_seqs.create())
		self.lines.reparentTo(self.render)
		
	def _voro_draw_alpha(self, alpha):
		new_line_seqs_thin = LineSegs()
		new_line_seqs_thin.setThickness(1.0) #TODO
		new_line_seqs_thin.setColor(Vec4(1,1,0,0.5)) #TODO			
		new_line_seqs_thick = LineSegs()
		new_line_seqs_thick.setThickness(2.0) #TODO
		new_line_seqs_thick.setColor(Vec4(0,0.9,0,0.75)) #TODO						
								
		point_cloud = PointCloud()
		cell_polygon = ConvexPolygons()
		make_face = False	
		n = 0
		for (points, faces) in self.input.current.voro:		
			if self.selected_cells[n]:				
				polygon_base_id = cell_polygon.vertex_count
				for (x,y,z,point_data) in points:	
					cell_polygon.add_vertex((x, y, z), 0.5, (1.0, 0.5, 0.2));
					if (point_data['alpha'] <= alpha):					
						point_cloud.add_vertex((x, y, z), 10.0, (0.1, 1, 0.1, 0.75))			
				
				for (poly, face_data) in faces:
					if (face_data['alpha'] <= alpha):
						make_face = True
					else:
						make_face = False
				
					(x,y,z,point_data) = points[poly[0][0]]					
					if (poly[0][1]['alpha'] > alpha):
						current_line_seq = new_line_seqs_thin
					else:
						current_line_seq = new_line_seqs_thick;
					current_line_seq.moveTo(Point3(x, y, z))
					
					if (make_face):
						cell_polygon.start_side()
						cell_polygon.add_side_point(polygon_base_id + poly[0][0])
						
					for i in xrange(1,len(poly)):									
						(x,y,z,point_data) = points[poly[i][0]]
						if (make_face):
							cell_polygon.add_side_point(polygon_base_id + poly[i][0])
						current_line_seq.drawTo(Point3(x, y, z))
						if (poly[i][1]['alpha'] > alpha):
							current_line_seq = new_line_seqs_thin
						else:
							current_line_seq = new_line_seqs_thick;
						current_line_seq.moveTo(Point3(x, y, z))
					(x,y,z,data) = points[poly[0][0]]
					current_line_seq.drawTo(Point3(x, y, z))
					
				if (self.grow_particles):
					self.balls[n].setColor(1.0, 1.0, 1.0, 0.25)	
				else:
					self.balls[n].setColor(1.0, 1.0, 1.0, 1.0)	
			else:
				if (self.grow_particles):
					self.balls[n].setColor(1.0, 1.0, 1.0, 0.05)	
				else:
					self.balls[n].setColor(1.0, 1.0, 1.0, 0.5)	
			n = n + 1

		if not(self.sides is None):
			self.sides.removeNode()
		#self.sides = self.render.attachNewNode(cell_polygon.create())
		#self.sides.reparentTo(self.render)	
		self.sides = cell_polygon.attach_to(self)
		self.sides.setTransparency(TransparencyAttrib.MDual, 1)		
		
		if not(self.lines_thin is None):
			self.lines_thin.removeNode()
		self.lines_thin = NodePath(new_line_seqs_thin.create())
		self.lines_thin.reparentTo(self.render)
		
		if not(self.lines_thick is None):
			self.lines_thick.removeNode()
		self.lines_thick = NodePath(new_line_seqs_thick.create())
		self.lines_thick.reparentTo(self.render)
		
		if not(self.vertices is None):
			self.vertices.removeNode()			
		#self.vertices = self.render.attachNewNode(point_cloud.create())
		#self.vertices.reparentTo(self.render)			
		self.vertices = point_cloud.attach_to(self)
		
		if not(self.lines is None):
			self.lines.removeNode()
			self.lines = None		


	def frame_iterator_task(self, task):			
		if self.need_update:
			balls_positions = self.input.current.balls;					
			for i, pos in enumerate(balls_positions):
				self.balls[i].setPos(pos.x, pos.y, pos.z)	
						
			if (self.grow_particles):
				for i, ball in enumerate(self.input.balls):
					r = sqrt(ball.r * ball.r + self.alpha)
					self.balls[i].setScale(r, r, r)
					self.balls[i].setTransparency(TransparencyAttrib.MDual, 1)		
					self.balls[i].setColor(1.0, 1.0, 1.0, 0.1)
			else:
				for i, ball in enumerate(self.input.balls):
					r = ball.r
					self.balls[i].setScale(r, r, r)
					self.balls[i].setTransparency(TransparencyAttrib.MDual, 1)		
					self.balls[i].setColor(1.0, 1.0, 1.0, 0.95)					
				
					
			if (self.show_shape):
				self._voro_draw_alpha(self.alpha)
			else:
				self._voro_draw()	
						
			self.need_update = False
		
		if self.is_running:		
			self.input.leap(1)		
			self.need_update = True
			#notify all obserwing functions :)
			for observer_func in self.observers:
				observer_func(self)
				
		return Task.cont

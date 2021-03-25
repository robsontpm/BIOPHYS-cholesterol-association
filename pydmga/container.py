import geometry
import dmga2py

#TODO: make get_particles(subset = None) for Container
#TODO: dodac klasy reprezentujace kinetic_structure

class ContainerIterator:
	'''
	used to iterate over Container like::
		
		for particle in container:
			pass
	
	:param container: Container for iteration
	'''
	def __init__(self, container):
		self._container = container
		self._container_handle = container._container_handle
		self._i = -1
		self._max = container.size()
		
	'''
	:return: self
	'''
	def __iter__(self):
		return self
	
	'''
	:return: next particle as a tuple of (id, x, y, z)
	'''
	def next(self):
		self._i += 1
		if (self._i >= self._max):
			raise StopIteration
		return dmga2py.basic_container_get(self._container_handle, self._i)

class Container:
	'''
	This holds particles in simulations
	'''
	def __init__(self, geometry):
		'''
		creates new container
		
		:param geometry: Geometry object to define the simulation box
		'''
		self.geometry = geometry;
		self._container_handle = dmga2py.new_basic_container(geometry._geometry_handle)
		self._size = 0
		
	def size(self):
		'''
		:return: number of particles in this container
		'''
		#TODO: maybe use the api?
		return self._size;
	
	#returns number of affected items
	def add(self, *args):
		'''
		Adds elements to this container.
		
		May work on many different elements:
		* add(tuple): tuple must contain 5 elements: (id, x, y, z, r), id is integer, others are double
		* add(iterable): iterable must contain tuples (as above)
		* add(id, z, y, z, r): elements as in tuple above
		* add(id, z, y, z): as above, but r is assumed to be 1.0
		
		Rises TypeError when elements are not tuples. 
		'''
		if len(args) == 1 and isinstance(args[0], tuple):
			(id, x, y, z, r) = args[0]
		elif len(args) == 1:
			try:
				for (id, x, y, z, r) in args[0]:
					self.add_by_coords(id, x, y, z, r)	
				return
			except TypeError:
				raise TypeError("Container::add: Not an iterable of tuples!")			
		elif len(args) == 5:
			id = args[0]
			x = args[1]
			y = args[2]
			z = args[3]
			r = args[4]
		elif len(args) == 4:		
			id = args[0]
			x = args[1]
			y = args[2]
			z = args[3]
			r = 1.0 #default
		else:
			raise TypeError("Container::add: cannot add that type of item!")		
		self.add_by_coords(id, x, y, z, r)
		
	def add_by_coords(self, id, x, y, z, r):
		'''
		for internal use, adds one element
		'''
		self._size += 1 #TODO: maybe use the api?
		(id, x, y, z, r) = self.geometry.transform(id, x, y, z, r)
		dmga2py.basic_container_add(self._container_handle, id, x, y, z, r)
		
	def add_raw(self, id, x, y, z, r):
		'''
		This adds one element without doing transformation on coordinates
		for boosting inserts, when we know that no transform is needed
		and we really need the speed... 
		'''	
		self._size += 1 #TODO: maybe use the api?
		dmga2py.basic_container_add(self._container_handle, id, x, y, z, r)		
		
		
	def get(self, number):
		'''
		Get the i-th inserted particle as a tuple (id, x, y, z)
		
		:param number: integer in range [0, size())
		:return: a particle with a given number
				
		**NOTE** number is not an ID! To get element by ID use find(ID)		
		'''
		return dmga2py.basic_container_get(self._container_handle, number)
		
	def find(self, id):
		'''
		Returns element with given ID
		
		:param id: is integer
		'''
		return dmga2py.basic_container_find(self._container_handle, id)
	
	def __del__(self):
		'''
		free all C++ objects if necessary
		'''	
		# print "freeing container_handle", self._container_handle
		dmga2py.free_object(self._container_handle) 
		
	def __iter__(self):
		'''
		makes this class iterable
		'''
		return ContainerIterator(self)
		

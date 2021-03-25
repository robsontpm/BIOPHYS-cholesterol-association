import dmga2py
from model import Cell  

#TODO: dodac klasy reprezentujace kinetic_structure

##################################################################################
##################################################################################
##################################################################################
class DiagramIterator:
	'''
	To iterate over cells in diagram.	
	Used internally in Diagram.__iter__()
	'''
	def __init__(self, diagram):
		'''
		create a new instance
		
		:param diagram: of class Diagram
		'''
		self._diagram = diagram
		self._i = -1;
		self._max = self._diagram.size()
	def __iter__(self):
		'''
		return self
		'''
		return self
	def next(self):
		'''
		returns next cell or rise StopIteration if no cells
		'''
		self._i += 1
		if (self._i >= self._max): 
			raise StopIteration
		return self._diagram.get_cell(self._i)

##################################################################################
##################################################################################
##################################################################################
class Diagram:
	'''
	This represents Voronoi Diagram of spheres (Power Diagram).
	You can iterate over all cells in diagram by::
	
		for cell in diagram:
			pass
			
	where cell is of class model.Cell
	'''
	def __init__(self, container, cache_on = False):
		'''
		Creates a new instance.
		
		:param containr: of class Container
		:param cache_on: turns on/off caching of computed cells. Default: False
		
		**NOTE** if cache is off (cache_on = False) in C++ version of the library user
		is responsible for freeing objects. That is if you need more memory, you can
		execute dmga2py.free_object on given particle (using handle). TODO: make it
		less cumbersome.
		'''
		self.container = container
		self._diagram_handle = dmga2py.new_diagram(container._container_handle, cache_on)
		
	def __iter__(self):
		'''
		to allow iteration of the form::
		
			for cell in diagram:
				pass
		'''
		return DiagramIterator(self)
	
	def size(self):
		'''
		:return: number of particles (cells) in this diagram
		'''
		return self.container.size()
		
	def get_cell(self, number):
		'''
		returns i-th cell (if necessary, compute it, or use cache)
		
		:param i: integer in range [0, size())
		:return: i-th cell as a object of class Cell		
		'''
		return Cell(dmga2py.diagram_get_cell(self._diagram_handle, number))  
		
	def get_cells(self, subset = None):
		'''
		Returns cells given by subset.
		
		:param subset: iterable that holds integers in range [0, size()) or None (default: None)
		
		#. **NOTE** if subset=None is given then all cells are returned.
		#. **NOTE** get_cells returns a list (not generator, may be heavy and time consuming!)
		#. **NOTE** for generating cells use iterator (for c in diagram: pass) 
		'''
		if (subset == None):
			return self.get_all_cells()
		result = []
		for i in subset:
			result.append(self.get_cell(i))		
		return result
		
	def get_all_cells(self):
		'''
		:return: all cells as list (may be heavy and time consuming!)
		
		**NOTE** for generating cells use iterator (for c in diagram: pass)
		'''
		return [cell for cell in self] # maybe you can make this faster by using C++ interface?
	
	def __del__(self):
		'''
		free all C++ objects if necessary
		'''	
		dmga2py.free_object(self._diagram_handle)
		

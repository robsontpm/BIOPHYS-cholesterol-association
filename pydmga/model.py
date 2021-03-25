import dmga2py
from math import sqrt

#TODO: implement Edges Collection and Iterator
#TODO: dodac klasy reprezentujace kinetic_structure


class EmptyIterator:
    '''
    Simple iterator that is allways empty 
    useful in some applications
    '''
    def __init__(self, *args):
        '''
        Ignores *args not needed
        '''
        pass

    def __iter__(self):
        '''
        standard
        '''
        return self

    def next(self):
        '''
        Always rises StopIteration
        '''
        raise StopIteration

##################################################################################
##################################################################################
##################################################################################
class CellVertex:
    '''
    This is a class that repreents single vertex of a Voronoi structure, that is
    the single point in space with some id, that can be used to create connection 
    relations (edges).
    
    :param _parent_cell_handle: handle to a C++ structure of a parent cell
    :param id: the number of this vertex (id) within a cell
    :param coords: tuple of coordinates (x,y,z)
    
    **NOTE** this class should not be created directly, you should access vertices by 
    routines in some Diagram.
    '''
    def __init__(self, _parent_cell_handle, id, coords):
        '''
        _parent_cell_handle is Handle
        id is int
        coords is tuple (double, double, double)
        '''
        self.id = id
        self._cell_handle = _parent_cell_handle
        (self.x, self.y, self.z) = coords
        
    def as_tuple(self):
        '''
        :return: this object as a tuple (id, x, y, z)
        '''
        return (self.id, self.x, self.y, self.z)
    
    def as_coords(self):
        '''
        :return: this object as a tuple (x, y, z)
        '''
        return (self.x,self.y,self.z)
    
    def __getitem__(self, key):
        '''
        allows to treat vertex as a tuple, i.e.
        to write v[0], v[1], v[2] for x, y, z coordinates respectively
        
        :param key: 0, 1, or 2, other -> KeyError
        :return: x, y or z respectively        
        '''
        if key == 0 or key == 'x': return x
        if key == 1 or key == 'y': return y
        if key == 2 or key == 'z': return z
        raise KeyError 

class CellVertexIterator:
    '''
    This is for internal use to allow iteration over Vertices in CellVertexCollection
    
    :param vertices: it is a CellVertexCollection
    :param as_coords: if True then the iteration will be returning tuples (x,y,z) of vertex rather than more involved CellVertex objects, default: False
    
    **NOTE** using ac_coords=True will be usually faster
    
    **NOTE** this should not be created directly, use CellVertexCollection from Cell
    to get this iterator.
    '''
    def __init__(self, vertices, as_coords = False):
        self._cell_handle = vertices._cell_handle
        self._max = vertices.size()
        self._i = -1;
        self._as_coords = as_coords

    def __iter__(self):
        '''
        returns self
        '''
        return self
    
    def next(self):
        '''
        returns next element
        '''
        self._i += 1
        if (self._i >= self._max): 
            raise StopIteration
        if (self._as_coords):
            return dmga2py.cell_get_vertex_coords(self._cell_handle, self._i);
        else:
            return CellVertex(self._cell_handle, self._i, dmga2py.cell_get_vertex_coords(self._cell_handle, self._i));

class CellVertexCollection:
    '''
    This represents iterable collection of vertices in a Voronoi cell
    
    this is usually accessible by Cell.vertices property.
    
    :param cell_handle: handle to C++ Cell counterpart
    
    **NOTE** this should not be created directly
    '''
    def __init__(self, cell_handle):
        self._cell_handle = cell_handle
        
    def size(self):
        '''
        :return: number of vertices in this cell
        '''
        return dmga2py.cell_get_vertex_count(self._cell_handle)
    
    def __iter__(self):
        '''
        allows iteration::
        
            for v in vertex_collection:
                pass
            
        where v is of class CellVertex
        '''
        try:     
            return CellVertexIterator(self)
        except:
            return EmptyIterator(self)
    
    def as_coords(self):
        '''
        allows iteration::
        
            for v in vertex_collection.as_coords():
                pass
                
        where v is simple tuple of coordinates (x,y,z)
        '''
        try:
            return CellVertexIterator(self, True)
        except:
            return EmptyIterator(self)

##################################################################################
##################################################################################
##################################################################################
class CellEdge:
    '''
    This class represents Edge of the Power Diagram
    :param parent_cell_handle: handle to C++ cell structure that contains this edge 
    :param id: unique id of this particular edge (each edge has unique id from 0 to cell.edges.size()-1)
    :param current: current = (v,j), v is first vertex index (ID) inside the cell, j is edge index at v
    :param inverse: current = (u,i), u is second vertex index (ID) inside the cell, i is edge index at u
    
    In inverse, i is the edge pointing backwards, i.e. from u to v
    
    We hold inverse to allow faster retrieval of inverse edge (useful)
    '''
    def __init__(self, parent_cell_handle, id, current, inverse):
        '''
        cell is Cell
        v is int - vertex index
        j is int - edge index at v
        u is int - next vertex (second vertex on edge j from v)
        i is int - edge from u to v at u (reverse of j)
        '''
        self.id = id
        self._cell_handle = parent_cell_handle
        (self._v, self._j) = current;
        (self._u, self._i) = inverse;
    
    def first(self):
        '''
        :return: v - id of the first vertex
        '''
        return self._v
    
    def second(self):
        '''
        :return: u - id of the second vertex
        '''
        return self._u
    
    def first_coords(self):
        '''
        :return: coordinates of the first vertex as a tuple (x,y,z)
        '''
        return dmga2py.cell_get_vertex_coords(self._cell_handle, self._v)
    
    def second_coords(self):
        '''
        :return: coordinates of the second vertex as a tuple (x,y,z)
        '''
        return dmga2py.cell_get_vertex_coords(self._cell_handle, self._u) 
    
    def as_pair_coords(self):
        '''
        :return: tuple of tuples: ((x_v, y_v, z_v), (x_u, y_u, z_u))
        
        may be used as a simple representation of an interval v -> u in 3D space
        
        :deprecated: this is replaced by as_coords  
        '''
        return (self.first_coords(), self.second_coords())
    
    def as_coords(self): #naming conventions...
        '''
        :return: tuple of tuples: ((x_v, y_v, z_v), (x_u, y_u, z_u))
        
        may be used as a simple representation of an interval v -> u in 3D space
        '''
        return (self.first_coords(), self.second_coords())
    
    def as_pair(self):
        '''
        :return: (v, u)
        '''
        return (self.first(), self.second())
    
    def as_tuple(self):
        '''
        :return: tuple (id, v, j, u, i)
        '''
        return (self.id, self._v, self._j, self._u, self._i)
    
class CellEdgesIterator:
    '''
    This is for internal use to allow iteration over Edges in CellEdgesCollection
    
    :param edges: reference to edges collection of class CellEdgeCollection 
    :param as_coords: if true, then next() will return pairs of coordinate tuples instead of CellEdge class instances, default: False
    
    **NOTE** this should not be created directly. Use cell.edges to get this iterator
    '''    
    def __init__(self, edges, as_coords = False):
        '''
        edges is CellEdgesCollection
        '''
        self._cell_handle = edges._cell_handle
        self._current = 0
        self._iter = dmga2py.cell_iterator(self._cell_handle)
        self._as_coords = as_coords
        
    def __iter__(self):
        '''
        returns self
        '''
        return self
    
    def next(self):
        '''
        go to the next element or throw StopIteration
        
        if as_coords = False then CellEdge class instance is returned, 
        otherwise tuple of coordinate tuples is returned.
        '''   
        if(dmga2py.cell_iterator_is_finished(self._iter) == 1):
            raise StopIteration
        
        dmga2py.cell_iterator_mark(self._iter)
        dmga2py.cell_iterator_forward(self._iter)
        if (dmga2py.cell_iterator_is_marked(self._iter) == 1):
            dmga2py.cell_iterator_jump(self._iter)
        if(dmga2py.cell_iterator_is_finished(self._iter) == 1):
            raise StopIteration
        current = self._current    
        self._current += 1
        dmga2py.cell_iterator_current_edge(self._iter)
        edge = CellEdge(self._cell_handle, current, dmga2py.cell_iterator_current(self._iter), dmga2py.cell_iterator_inverse(self._iter))
        if self._as_coords:
            return edge.as_coords()
        else:
            return edge 

class CellEdgesCollection:
    '''
    This represents iterable collection of edges in a Voronoi cell
    
    :param cell_handle: handle to associated C++ cell structure
    
    **NOTE** Should not be created directly. Use cell.edges to get this. 
    '''
    def __init__(self, cell_handle):
        '''
        cell is Cell
        '''
        self._cell_handle = cell_handle
        
    def size(self):
        '''
        :return: number of edges in this cell
        '''
        return dmga2py.cell_get_edges_count(self._cell_handle)
    
    def __iter__(self):
        '''
        allows iteration::
        
            for edge in edges_collection:
                pass
                
        where edge is of class CellEdge
        '''
        try:
            return CellEdgesIterator(self) 
        except:
            return EmptyIterator(self)
    
    def as_coords(self):
        '''
        allows iteration::
        
            for (vertex_1_coords, vertex_2_coords) in edges_collection.as_coords():
                pass
                
        where vertex_1_coords and vertex_2_coords are tuples (x,y,z)
        '''
        try:
            return CellEdgesIterator(self, True)
        except:
            return EmptyIterator(self)
    
##################################################################################
##################################################################################
##################################################################################
class CellSideIterator:
    '''
    This is for internal use to allow iteration over vertex coordinates in CellSide
    '''    
    def __init__(self, cell_side):
        '''
        cell_side is CellSide
        '''
        self._cell_side = cell_side
        self._vertices_iterator = cell_side.vertices.__iter__()
        
    def __iter__(self):
        '''
        returns self
        '''
        return self
    
    def next(self):
        '''
        returns coordinates of the next vertex as a tuple (x,y,z)
        '''
        index = self._vertices_iterator.next()
        return dmga2py.cell_get_vertex_coords(self._cell_side._cell_handle, index)
        

class CellSide:
    '''
    represents a single side of a Voronoi cell
    
    :param cell_handle: handle to associated C++ cell structure
    :param number: the number of the face (unique ID from 0 to cell.sides.size()-1)
    :param neighbour: number of the neighbour who created this side
    :param vertices: list of vertices on this side
    
    
    **NOTE** Should not be created directly. Use cell.sides to get this.     
    '''
    def __init__(self, cell_handle, number, neighbour, vertices):
        '''
        creates new side
        :param cell_handle: handle to associated C++ cell structure
        :param number: the number of the face (unique ID from 0 to cell.sides.size()-1)
        :param neighbour: number of the neighbour who created this side
        :param vertices: list of vertices on this side
        '''
        self._cell_handle = cell_handle
        self.number = number
        self.neighbour = neighbour
        self.vertices = vertices
        
    def size(self):
        '''
        returns number of vertices on this side
        '''
        return len(self.vertices)
        
    def __iter__(self):
        '''
        allows iteration over this side vertices::
        
            for (x,y,z) in cell.side:
                pass
                
        TODO: should it return self.vertices.__iter__() instead?
        '''
        try:
            return CellSideIterator(self)
        except:
            return EmptyIterator(self)
    
    def as_coords(self):
        '''
        :return: this side as list of coordinates (x,y,z) of all vertices
        '''
        return [c for c in self.__iter__()]
    
    def as_list(self):
        '''
        :return: this side as list of vertex indices
        '''
        return self.vertices
    
    def __del__(self):
        '''
        free all C++ objects if necessary
        ''' 
        pass
    
    def area(self):
        '''
        :return: area of this side
        '''
        area = 0.0
        iter = self.__iter__()
        v1 = iter.next() 
        try:            
            v2 = iter.next()
            v3 = iter.next()
            while(True):
                u1 = [v2[i]-v1[i] for i in range(3)]
                u2 = [v3[i]-v1[i] for i in range(3)]
                wx = u1[1]*u2[2] - u1[2]*u2[1]; # cross product
                wy = u1[2]*u2[0] - u1[0]*u2[2]; # cross product
                wz = u1[0]*u2[1] - u1[1]*u2[0]; # cross product
                area += sqrt(wx**2 + wy**2 + wz**2);
                v2 = v3
                v3 = iter.next()
        except StopIteration:
            pass
        if (area < 0): area = -area
        return area / 2.0
        
class CellSidesIterator:
    '''
    This is for internal use to allow iteration over Sides in CellEdgesCollection
    it returns sides as lists of CellSide objects or as lists of coordinates.
    
    :param sides: is CellSidesCollection
    :param as_coords: is Boolean and controls if iterator should return coords (True) or CellSide class instances (False) (default: False)
    
    '''    
    def __init__(self, sides, as_coords = False):
        '''
        :param sides: is CellSidesCollection
        :param as_coords: is Boolean and controls if iterator should return coords (True) or CellSide class instances (False) (default: False)
        it is internally used in CellSidesCollection.__iter__() to allow iteration
        '''
        self._current = 0
        self._cell_handle = sides._cell_handle
        self._iter = dmga2py.cell_iterator(self._cell_handle)
        self._as_coords = as_coords
        
    def __iter__(self):
        '''
        standard
        '''
        return self
    
    def next(self):
        '''
        go to the next element or throw StopIteration
        '''        
        if(dmga2py.cell_iterator_is_finished(self._iter) == 1):
            raise StopIteration
        side = []
        # TODO: rename cell_iterator_current_side_id to cell_iterator_current_side_neighbour
        neighbour = dmga2py.cell_iterator_current_side_neighbour_id(self._cell_handle, self._iter)
        while (dmga2py.cell_iterator_is_finished(self._iter) == 0):
            side.append(dmga2py.cell_iterator_current_vertex(self._iter))
            dmga2py.cell_iterator_mark(self._iter)
            dmga2py.cell_iterator_forward(self._iter)
            if (dmga2py.cell_iterator_is_marked(self._iter) == 1):
                dmga2py.cell_iterator_jump(self._iter)
                break
        current = self._current
        self._current += 1
        cell_side = CellSide(self._cell_handle, current, neighbour, side)
        if (self._as_coords):
            return cell_side.as_coords()
        else:
            return cell_side
        

class CellSidesCollection:
    '''
    This represents iterable collection of sides in a Voronoi cell
    '''    
    def __init__(self, cell_handle):
        self._cell_handle = cell_handle       
        
    def size(self):
        '''
        returns number of sides (neighbours) in this cell
        '''
        return dmga2py.cell_get_sides_count(self._cell_handle)
    
    def __iter__(self):
        '''
        allows iteration::
        
            for side in side_collection:
                pass
                
        where side is of class CellSide
        '''
        try:
            return CellSidesIterator(self)  
        except:
            return EmptyIterator(self)
    
    def as_coords(self):
        '''
        allows iteration::
        
            for side in side_collection:
                pass
                
        where side is a list of coordinate tuples (x,y,z) of consecutive vertices.
        '''
        try:
            return CellSidesIterator(self, True)
        except:
            return EmptyIterator(self)


##################################################################################
##################################################################################
##################################################################################
class Cell:
    '''
    Represents Voronoi Cell
    
    It is created by Diagram class
    
    :param cell_handle: handle to C++ object representing the cell
    
    it has three public attributes, which are collections of sides, edges and vertices. 
    
    It allows for constructs such as::
    
        for side in cell.sides:
            print(side.area())
            print(side.as_coords())
        
        print(cell.sides.size())
        
        for edge in cell.edges:
            print(edge.as_coords())
            print(edge.id)
        
        print cell.edges.size()
        
        for vertex in cell.vertices:
            print(vertex.as_coords())
            print(vertex.id)
        
        print(cell.vertices.size())
                
    '''
    def __init__(self, cell_handle):
        '''
        for internal use, it accepts cell_handle to C++ Cell class object
        '''
        self._cell_handle = cell_handle
        self.vertices = CellVertexCollection(self._cell_handle)
        self.edges = CellEdgesCollection(self._cell_handle)
        self.sides = CellSidesCollection(self._cell_handle)

    def is_empty(self):
        return True if dmga2py.cell_is_empty(self._cell_handle) else False
        
    def volume(self):
        '''
        :return: volume of this cell
        '''
        return dmga2py.cell_get_volume(self._cell_handle)
    
    def area(self):
        '''
        :return: area of this cell
        '''
        return dmga2py.cell_get_area(self._cell_handle)
    
    def __iter__(self):
        '''
        allows iteration on sides
        '''
        if self.is_empty():
            return EmptyIterator(self)
        return self.sides.__iter__()   
       
    def get_sides_as_coords(self):
        '''
        returns all sides as one list, thus may be time and space consuming
        the resulting list contains only lists of 3D coordinates of vertices
        of each side.
        for generating better use sides iterable collection::
        
            for side in cell.sides:
                pass
                 
        '''
        if self.is_empty():
            []
        result = []
        iter = dmga2py.cell_iterator(self._cell_handle)
        side = []
        while (dmga2py.cell_iterator_is_finished(iter) == 0):
            side.append(dmga2py.cell_iterator_current_vertex_coords(iter))
            dmga2py.cell_iterator_mark(iter)
            dmga2py.cell_iterator_forward(iter)
            if (dmga2py.cell_iterator_is_marked(iter) == 1):
                result.append(side);
                side = []
                dmga2py.cell_iterator_jump(iter)
        #TODO: przemyslec, czy mam usuwac iterator handle czy nie!  
        return result
    
    def get_sides(self):
        '''
        returns all sides as one list, thus may be time and space consuming
        the resulting list contains CellSide objects.
        for generating better use sides iterable collection::
        
            for (neighbour, side) in cell.sides:
                pass
                
        where neighbour is ??? and side is ???                
        '''    
        if self.is_empty():
            return []   
        result = []
        iter = dmga2py.cell_iterator(self._cell_handle)
        side = []
        current = 0
        while (dmga2py.cell_iterator_is_finished(iter) == 0):
            side.append(dmga2py.cell_iterator_current_vertex_coords(iter))
            neighbour = dmga2py.cell_iterator_current_side_id(self._cell_handle, iter)
            dmga2py.cell_iterator_mark(iter)
            dmga2py.cell_iterator_forward(iter)
            if (dmga2py.cell_iterator_is_marked(iter) == 1):
                result.append( CellSide(self, current, neighbour, side) );
                current += 1
                side = []
                dmga2py.cell_iterator_jump(iter)
        #TODO: przemyslec, czy mam usuwac iterator handle czy nie!
        return result
    
    def __del__(self):
        '''
        free all C++ objects if necessary
        '''    
        # print "freeing cell_handle", self._cell_handle
        dmga2py.free_object(self._cell_handle) 
        
   

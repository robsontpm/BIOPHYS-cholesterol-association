import dmga2py
from diagram import Diagram
from container import Container
from model import Cell, CellVertex, CellEdge, CellSide    

#TODO: dodac klasy reprezentujace kinetic_structure

class ShapeIterator:
    '''
    To iterate over cells in a Shape.
    Used internally in Shape.__iter__()
    
    **NOTE** This is default Shape iterator and should work for
    all the shapes. It is necessary for the Shape to have get_cell(i) function
    '''
    def __init__(self, shape):
        '''
        create a new instance
        '''
        self._shape = shape
        self._i = -1;
        self._max = self._shape.size()
    def __iter__(self):
        '''
        return self
        '''
        return self
    def next(self):
        '''
        returns next cell shape or rise StopIteration if no cells
        '''
        self._i += 1
        if (self._i >= self._max): 
            raise StopIteration
        return self._shape.get_cell(self._i)

class Shape:  
    '''
    This is abstract Shape class, it is used to define various properties of
    the voronoi diagram that can be inferred using distance and neighbouring information.
    Derived class should implement size() and get_cell(i) functions.
    
    Shape provides standard iterator and get_cells(subset) functions.
    
    allows iteration in the form::
        
        for shape_cell in shape:
            pass
            
    **NOTE** This is abstract class and should not be instantiated
    '''  
    def __iter__(self):
        '''
        allows iteration in the form::
        
            for shape_cell in shape:
                pass
                
        where shape_cell is of class AlphaShapeCell
        '''
        return ShapeIterator(self)
    
    def size(self):
        '''
        Not implemented yet
        '''
        raise NotImplementedError("This is abstract class!")
    
    def get_cell(self, number):
        '''
        Not implemented yet
        '''
        raise NotImplementedError("This is abstract class!")
    
    def get_cells(self, subset = None):
        '''
        returns list of cells for given subset
        subset is iterable of integers or None (default: None)
        if None is given then all cells will be returned
        this returns list so may be time and space consuming
        for generating use iteration::
        
            for shape_cell in shape:
                pass
            
        where shape_cell is of class AlphaShapeCell
        '''
        if subset is None:
            return [shape_cell for shape_cell in self]
        else:
            result = []
            for i in subset:
                results.append(slef.get_cell(i))
            return result 

##################################################################################
# Alpha Shapes ###################################################################
##################################################################################
class AlphaSimplex:
    '''
    represents simple element of Power Diagram with given alpha in the coressponding alpha
    complex. As we construct Voronoi Diagram, we do not have literally the Alpha Complex (which is
    subset of Delaunay Triangulation) but we can use duality. Of course we have duality only in case of
    General Position of balls, but for us we only care for the time when given element of Power Diagram is
    touched by a growing ball, so we can treat equally elements for balls not in general position as in 
    the case of standard elements.
    
    We also need to note the coresspondence (PG is Power Diagram, RT is regular triangulation (Delaunay)):
    * VERTEX in PG <-> Tetrahedron in RT
    * EDGE in PG <-> Triangle in RT
    * SIDE in PG <-> Edge in RT
    
    :param alpha_shape_cell_handle: is a handle to DMGA C++ object which holds the AlphaShapeCell structure
    :param dim: is the dimension of this simplex
    :param v: is a base vertex for the this simplex. For vertex it is simply this vertex, for Edge it is base vertex, for side it is starting vertex.
    :param j: is edge number at v at which this simplex starts, for vertex it should be ignored, for edge it is the edge, for side it is the first edge on this side from v
    :param alpha_threshold: is alpha value at which this simplex is touched by the ball with radius^2 = r^2 + alpha_threshold
    :param id: is the inner id within the cell of this simplex, for vertex it is this vertex number, for edge it is id of the edge (unique) and for side it is the number of side with standard iteration
    :param is_killed: is a flag that tells if this alpha_threshold is not the smallest power distance form ball center, but other alpha of other simplex. It happens when nearest point doesn't belong to diagram part of this simplex.
    
    * **NOTE** In the structure j will be stored in edge property.    
    * **NOTE** In the structure v will be stored in vertex property.    
    * **NOTE** In the structure dim will be stored in dimension property.  
    
    **NOTE** This class should not be instantiated directly - you should use AlphaShapeCell iteration ability or accessors to get it  
    '''
    DIM_VERTEX = 3  #3-d simplex is vertex as dually in Delaynay triangulation it vould be Tetrahedron
    DIM_EDGE = 2    #2-d simplex is edge as dually in Delaynay triangulation it vould be Triangle (side)
    DIM_SIDE = 1    #1-d simplex is edge as dually in Delaynay triangulation it vould be Edge connecting two particles (nearest neighbours)
    #we do not use 0-d simplexes as those in Delaunay simplex are vertices - particles defining this triangulation and they are allwejs in Alpha-Complex
    
    def __init__(self, alpha_shape_cell_handle, dim, v, j, alpha_threshold, id, is_killed):
        '''
        dimension (dim) is the dimension of this simplex
        vertex is a base vertex for the this simplex. For vertex it is simply this vertex, for Edge it is base vertex, for side it is starting vertex.
        edge is edge number at v at which this simplex starts, for vertex it should be ignored, for edge it is the edge, for side it is the first edge on this side from v
        alpha_threshold is alpha value at which this simplex is touched by the ball with radius^2 = r^2 + alpha_threshold
        id is the inner id within the cell of this simplex, for vertex it is this vertex number, for edge it is id of the edge (unique) and for side it is the number of side with standard iteration
        is_killed is a flag that tells if this alpha_threshold is not the smallest power distance form ball center, but other alpha of other simplex. It happens when nearest point doesn't belong to diagram part of this simplex.  
        '''
        self._alpha_shape_cell_handle = alpha_shape_cell_handle
        self.dimension = dim
        self.vertex = v
        self.edge = j
        self.alpha_threshold = alpha_threshold
        self.id = id
        self.is_killed = is_killed
        
    def diagram_part(self):
        '''
        :return: corresponding part of the diagram
        
        Returned value is either vertex as model.CellVertex, side as model.CellSide, or edge model.CellEdge 
        '''
        (cell_handle, complex_handle) = dmga2py.alpha_shape_cell_as_tuple(self._alpha_shape_cell_handle)
        (v, j) = (self.vertex, self.edge)
        if (self.dimension == self.DIM_VERTEX): 
            #we have vertex in Voronoi (tetrahedron in Delaunay)
            return CellVertex(cell_handle, self.id, dmga2py.cell_get_vertex_coords(cell_handle, v))
        if (self.dimension == self.DIM_EDGE):
            #we have edge in Voronoi (side in Delaunay)
            return CellEdge(cell_handle, self.id, (v, j), dmga2py.cell_edge_get_inverse(cell_handle, v, j))
        if (self.dimension == self.DIM_SIDE): 
            #we have side in Voronoi (edge i Delaunay)
            side_neighbour = dmga2py.cell_edge_get_neighbour_id(cell_handle, v, j)
            (u, i) = dmga2py.cell_edge_get_next(cell_handle, v, j)
            side_id = self.id
            side = [v]
            while (u, i) != (v, j): # this suffice, as we are interested only on this side
                side.append(u)
                (u, i) = dmga2py.cell_edge_get_next(cell_handle, u, i)
            return CellSide(cell_handle, side_id, side_neighbour, side)
        raise ValueError("AlphaSimplex.diagram_part(): simplex with dimension " + str(self.dimension))
    
    def as_tuple(self):
        '''
        returns this element as a tuple::
        
            (dimension, vertex, edge, alpha_threshold, id, is_killed)
        '''
        return (self.dimension, self.vertex, self.edge, self.alpha_threshold, self.id, self.is_killed)
    
    def as_dict(self):
        '''
        returns this element as a dictionary::
        
            "dimension": self.dimension, 
            "vertex": self.vertex, 
            "edge": self.edge, 
            "alpha_threshold": self.alpha_threshold, 
            "id": self.id, 
            "is_killed": self.is_killed
            
        '''
        return {"dimension": self.dimension, 
                "vertex": self.vertex, 
                "edge": self.edge, 
                "alpha_threshold": self.alpha_threshold, 
                "id": self.id, 
                "is_killed": self.is_killed}
        
    def as_dict_short(self):
        '''
        returns this element as a dictionary::
        
            "d": self.dimension, 
            "v": self.vertex, 
            "e": self.edge, 
            "a": self.alpha_threshold, 
            "i": self.id, 
            "k": self.is_killed
            
        '''        
        return {"d": self.dimension, 
                "v": self.vertex, 
                "e": self.edge, 
                "a": self.alpha_threshold, 
                "i": self.id, 
                "k": self.is_killed}
        
class AlphaShapeCellIterator:
    '''
    To iterate over all simplexes in given cell
    used internally in AlphaShapeCell.shape()
    
    The iteration is done alpha-wise,m that is from the smallest alpha to the biggest.
    
    :param alpha: Float, this is the upper bund on the elements we want to iterate, so only elements with element.alpha <= alpha will be present.
    :param start: iteration starts here, it is sometimes useful. 
    
    example::
    
        AlphaShapeCellIterator(shape_cell, 5.0, 10)
        
    will create an iterator for all elements with alpha <= 5.0. Suppose there are 100 such elements. Then
    the iteration will start from the 10-th element. 
    
    **NOTE** this class should not be instantiated directly - use AlphaCellShape iteration to get it.
    
    '''
    def __init__(self, alpha_shape_cell, alpha, start = 0):
        '''
        create a new instance
        '''
        self._alpha = alpha
        self._alpha_shape_cell = alpha_shape_cell
        self._i = start-1
        self._max = dmga2py.alpha_shape_complex_upper_bound(self._alpha_shape_cell._complex_handle, alpha)
        
    def __iter__(self):
        '''
        return self
        '''
        return self
    
    def next(self):
        '''
        returns next cell simplex or rise StopIteration if no cells
        '''
        self._i += 1
        if self._i >= self._max:
            raise StopIteration
        (dim,v,j,alpha,id,killed) = dmga2py.alpha_shape_complex_get_simplex(self._alpha_shape_cell._complex_handle, self._i)
        return AlphaSimplex(self._alpha_shape_cell._shape_cell_handle, dim, v, j, alpha, id, killed)
  
class AlphaShapeCell(Cell):
    '''
    its the extension of a base Cell, it provides the same operations as Cell (for convenience)    
    and defines an alpha shape of a single cell, i.e. the sorted 
    list of 1, 2 or 3-D simplexes with their alhpa threshold - 
    value of alpha at which alpha ball with r'^2 = r^2 + alpha 'kills'
    this simplex (i.e. when intersection of ball and siplex is nonempty).
    In our computation we use dual of the Delaunay triangulation - Power diagram, but
    we also allow for non-general positions of balls.
    The shape is accessible by shape(alpha) method, if no alpha specified then we assume alpha = max_alpha_threshold
    
    **NOTE** This class should not be instantiated directly - it is returned from AlphaShape by access functions.  
    
    :param shape_cell_handle: tells which shape_cell C++ handle to use
    :param cell_handle: tells which Voronoi cell is used (for which Voronoi Cell), this is C++ handle
    :param complex_handle: tells which complex is used (should be relative to cell_handle), this is C++ handle
    ''' 
    def __init__(self, shape_cell_handle, cell_handle, complex_handle):
        Cell.__init__(self, cell_handle)
        self._shape_cell_handle = shape_cell_handle
        self._cell_handle = cell_handle
        self._complex_handle = complex_handle
    
    def shape(self, alpha = None):
        '''
        allows for iteration over all simplices in shape::
        
            for simplex in alpha_cell.shape():
                pass
                
        or only for those simplices that has alpha lower than given alpha threshold::
        
            for simplex in alpha_cell.shape(10.0):
                #iterates only over simplices with alpha <= 10.0
                pass
                
        '''
        if (alpha is None):
            alpha = self.max_alpha_threshold()
        return AlphaShapeCellIterator(self, alpha)
    
    def __del__(self):
        '''
        free all C++ objects if necessary
        '''    
        pass
    
    def max_alpha_threshold(self):
        '''
        :return: max alpha_threshold - value of alpha at which all elements from this cell belongs to alpha complex
        '''        
        #return dmga2py.alpha_shape_cell_get_max_alpha_threshold(self._shape_cell_handle)
        dmga2py.alpha_shape_complex_print(self._complex_handle)    
        return dmga2py.alpha_shape_complex_get_max_alpha_threshold(self._complex_handle)
    
class AlphaShape(Shape):
    '''
    This class represents a description of the Alpha Shape of a given Diagram
    
    :param diagram_or_container: Diagram or Container class to use for computations
    
    **NOTE** if Container class is used to instantiate then the class will create Diagram for internal use
    '''
    def __init__(self, diagram_or_container):
        '''
        diagram_or_container is either Diagram or Container, if container, then Diagram is internally created
        '''
        # make sure we have diagram
        if isinstance(diagram_or_container, Diagram):
            # print "AlphaShape::__init__(): we have Diagram"
            self.diagram = diagram_or_container            
        elif isinstance(diagram_or_container, Container):
            # make new diagram
            # print "AlphaShape::__init__(): we have Container"
            self.diagram = Diagram(diagram_or_container)
        else:
            raise TypeError("AlphaShape.__init__(): Diagram or Container required but something else given")        
        # do rest of the initialization
        self._diagram_handle = self.diagram._diagram_handle
        # print "Diagram handle =", self.diagram._diagram_handle
        # create item on the C++ side
        self._shape_handle = dmga2py.new_alpha_shape(self._diagram_handle)
    
    def size(self):
        '''
        :return: number of cells in this shape, the same as size in diagram
        '''
        return self.diagram.size()
    
    def get_cell(self, number):
        '''
        :return: one cell of this shape
        
        The class of shape depends on the type of the shape
        in Alpha Shape it is Alpha Complex (subset of...) for this cell
        '''
        (shape_cell_handle, cell_handle, complex_handle) = dmga2py.alpha_shape_get_shape_cell(self._shape_handle, number)
        # complex_handle = dmga2py.alpha_shape_get_complex(self._shape_handle, number)
        return AlphaShapeCell(shape_cell_handle, cell_handle, complex_handle)
    
    def __del__(self):
        '''
        free all C++ objects if necessary
        '''    
        dmga2py.free_object(self._shape_handle)  
        
    def __iter__(self):
        '''
        :return: iterator over all shape cells in this shape
        '''
        return ShapeIterator(self)    
    
    def max_alpha_threshold(self):
        '''
        :return: max alpha_threshold - value of alpha at which all elements from diagram belongs to the alpha-complex for all cells
        
        **NOTE** max_threshold is computed only for those cells that were actually created by for 
        example get_cell() or get_cells() or if diagram has cache_on property set to ON (True).
        '''        
        return dmga2py.alpha_shape_get_max_alpha_threshold(self._shape_handle)         
        
##################################################################################
# SASA Shapes ####################################################################
##################################################################################
class SASAArc:
    '''
    This class describes one arc of the SASA shape, that is an arc on the 
    surface of the sphere. 
    
    :param first:first end point of the line
    :param on_plane: the central point of the arc lying on the Voronoi plane (plane containing some side) 
    :param on_sphere: the central point on the surface of the sphere
    :param second: second end point
    
    on_plane and on_sphere points may be used to calculate various quentities of the arc (angles, areas, etc)
    and to draw the arc. 
    '''
    def __init__(self, first, on_plane, on_sphere, second):
        self.first = first
        self.on_plane = on_plane
        self.on_sphere = on_sphere
        self.second = second
        
    def __str__(self):
        return "Arc<{0}, <{1}, {2}>, {3}>".format(self.first, self.on_plane, self.on_sphere, self.second)
    
class SASAContoursIterator:
    '''
    To iterate over all simplexes in given cell
    used internally in AlphaCell.__iter__()
    
    :param contour_iterator_handle: C++ handle to an iterator
    '''
    def __init__(self, contour_iterator_handle):
        '''
        create a new instance
        '''
        self._iter_handle = contour_iterator_handle #dmga2py.sasa_contours_new_border_iterator(sasa_cell._contour_handle)
        if (dmga2py.sasa_contours_iterator_ready(self._iter_handle)):
            self._stop = False
        else:
            self._stop = True
        #self._iter_handle = sasa_cell_iterator_handle
        #self._doing = dmga2py.sasa_shape_cell_iterator_ready(sasa_cell_iterator_handle) 
        
    def __iter__(self):
        '''
        returns self
        '''
        return self
    
    def next(self):
        '''
        returns next contour
        '''
        if (self._stop):
            raise StopIteration
        (first, on_plane, on_sphere, second) = dmga2py.sasa_contours_iterator_as_tuple_coords(self._iter_handle)
        if (dmga2py.sasa_contours_iterator_iterate(self._iter_handle) == -1):
            self._stop = True
        return SASAArc(first, on_plane, on_sphere, second)

class SASACell(Cell):
    '''
    This represents a SASA contours description of a given cell. Those contours may be used
    to compute excluded/included area and volume. The excluded area is an area of the sphere
    surface outside the voronoi cell, included is the inverse. The same for the volume. 
    Sum of all included area is a SAS Area of the system.  
    
    TODO: this should be SASShape? (naming) 
    '''
    def __init__(self, shape_cell_handle, cell_handle, contour_handle, number):
        Cell.__init__(self, cell_handle)
        self._sasa_cell_handle = shape_cell_handle
        self._cell_handle = cell_handle
        self._contour_handle = contour_handle        
    
    def sas_area(self):
        '''
        :return: included area - the area of the sphere surface contained inside the cell
        '''
        return dmga2py.sasa_shape_cell_get_area(self._sasa_cell_handle)
    
    def sas_volume(self):
        '''
        :return: included volume - the volume of the ball volume contained inside the cell
        '''
        dmga2py.sasa_shape_cell_get_volume(self._sasa_cell_handle)
        
    def excluded_area(self):
        '''
        :return: excluded area - the area of the sphere surface outside the cell
        '''
        return dmga2py.sasa_shape_cell_get_excluded_area(self._sasa_cell_handle)
    
    def excluded_volume(self):
        '''
        :return: excluded volume - the volume of the ball outside the cell
        '''
        dmga2py.sasa_shape_cell_get_excluded_volume(self._sasa_cell_handle) 
    
    def border_arcs(self):        
        '''
        return: iterator over all border arcs
        
        Border arcs are those arcs that form the border of the SAS Contour on the surface of 
        the sphere. It is a result of intersection of the Sphere with the Voronoi Cell.
        
        This iterator should be used when one needs to draw (visualize) the contour 
        '''
        return SASAContoursIterator(dmga2py.sasa_contours_new_border_iterator(self._contour_handle))
    
    def polygon_arcs(self):
        '''
        return: iterator over all polygon arcs
        
        Polygon arcs are those arcs that forms the polygons inside the excluded area. Notice that
        border arcs form closed curve consisting of parts of circles (thus name arcs was given). Joining 
        each point connecting two different arcs to the centers of the balls on the sphare surface with great 
        arcs we get parts of the domes and the rest is a (sum of) spherical polygons. 
        The border of spherical polygons consists of parts of those great arcs from the point to the 
        center on sphere. The center on sphere is a pint on sphere where the line connecting two neighbouring
        particles crosses the sphere.
        
        '''
        return SASAContoursIterator(dmga2py.sasa_contours_new_border_iterator(self._contour_handle))
    
    def shape(self, mode = 0):
        '''
        returns the shape iterator over this SASAShape. 
        :param mode: defines which shape to return, default is Border Arcs (mode=0), the other is polygon_arcs (mode=1) 
        '''
        if (mode == 0):
            return self.border_arcs()
        return self.polygon_arcs()
    
    def __iter__(self):
        return Cell.__iter__(self)

class SASAShape(Shape):
    '''
    This is a Solvent Accessible Surface Shape (SAS) for a given Diagram. It is
    used mainly to compute Solvent Accessible Surface Area (SASA)
    
    SAS is a shape of the union of balls and thus is strongly connected with Alpha SHape.
    Alpha SHape is more like a topological description of the Shape (what is connected with what), while
    SAS is the description of the actual union of balls. If you are interested only in topological aspect
    we suggest Alpha Shapes as those are a bit faster to compute.  
    
    :param diagram_or_container: Diagram or Container class to use for computations
    
    **NOTE** if Container class is used to instantiate then the class will create Diagram for internal use    
    '''
    def __init__(self, diagram_or_container):
        '''
        diagram_or_container is either Diagram or Container, if container, then Diagram is internally created
        '''
        # make sure we have diagram
        if isinstance(diagram_or_container, Diagram):
            self.diagram = diagram_or_container            
        elif isinstance(diagram_or_container, Container):
            # make new diagram
            self.diagram = Diagram(diagram_or_container)
        else:
            raise TypeError("AlphaShape.__init__(): Diagram or Container required but something else given")        
        # do rest of the initialization
        self._diagram_handle = self.diagram._diagram_handle
        # create item on the C++ side
        self._shape_handle = dmga2py.new_sasa_shape(self._diagram_handle)
    
    def size(self):
        '''
        :return: number of cells in this shape, the same as size in diagram
        '''
        return self.diagram.size()
    
    def get_cell(self, number):
        '''
        :return: one cell of this shape
        
        the class of shape depends on the type of the shape
        in SASA Shape it is SASACell for this cell that holds the COntour information (SASArcs)
        '''
        (shape_cell_handle, cell_handle, contour_handle) = dmga2py.sasa_shape_get_shape_cell(self._shape_handle, number)
        return SASACell(shape_cell_handle, cell_handle, contour_handle, number)
    
    def __del__(self):
        '''
        free all C++ objects if necessary
        ''' 
        # TODO: why it is caousing core dump??????? 
        dmga2py.free_object(self._shape_handle)
        pass
        
    def __iter__(self):
        '''
        :return: iterator over all shape cells in this shape
        '''
        return ShapeIterator(self)
    
    def sas_area(self):
        '''
        :return: area of all computed cells in this SASA Shape that is area of the intersection of sphere with its voronoi cell 
        '''
        return dmga2py.sasa_shape_get_area(self._shape_handle)
    
    def sas_volume(self):
        '''
        :return: volume of all computed cells in this SASA Shape that is volume of the intersection of ball with its voronoi cell 
        '''
        return dmga2py.sasa_shape_get_volume(self._shape_handle)    
    
    
    
     
import bempp
import numpy as np
from binary_search import coordBinarySearch

def merge (grid1, grid2, tolerance):
    """
    Merge two meshes. Found duplicate nodes. No creation of new elements.

    Parameters
    ----------
    grid1 : bempp.grid.Grid object
    First mesh
    grid2 : bempp.grid.Grid object
    Second mesh
    tolerance : float
    Tolerance for searching of duplicated nodes

    """

    # read vertices of two meshes as tuples
    vertex_iterator = grid1.leaf_view.entity_iterator(2)
    index_set = grid1.leaf_view.index_set()

    ord1 = tuple([index_set.entity_index(vertex),vertex.geometry.corners[:,0]]
        for vertex in vertex_iterator)       

    vertex_iterator = grid2.leaf_view.entity_iterator(2)
    index_set = grid2.leaf_view.index_set()

    ord2 = tuple([index_set.entity_index(vertex),vertex.geometry.corners[:,0]]
        for vertex in vertex_iterator)

    #number of nodes in two meshes
    number_of_nodes1 = len(ord1);
    number_of_nodes2 = len(ord2);

    # sort vertices in second mesh by x-coord
    sortArrX = tuple( sorted( ord2, key = lambda t: t[1][0] ))
    transposeX = zip(*sortArrX)
    coordX = zip(*(transposeX[1]))[0]

    # find duplicates
    duplicateIndices = [] #renumeration of duplicate nodes in mesh 2

    for i1 in xrange(0,number_of_nodes1):
        duplicateIndex = coordBinarySearch (sortArrX, coordX, ord1[i1][1], tolerance)

        if (duplicateIndex > -1):
            duplicateIndices.append([duplicateIndex,ord1[i1][0]])     

    n_duplicated = len(duplicateIndices)  
    print ("duplicated nodes: ",len(duplicateIndices))

    # merge

    #read vertices/elements from first mesh
    vertices_one = grid1.leaf_view.vertices;
    elements_one = grid1.leaf_view.elements; 
    domain_indices_one = grid1.leaf_view.domain_indices;

    #read vertices/elements from second mesh
    vertices_two = grid2.leaf_view.vertices;
    elements_two = grid2.leaf_view.elements;
    domain_indices_two = grid2.leaf_view.domain_indices;

    # merge nodes
    vertices_merged = vertices_one.tolist()
    n_renumerated = 0
    duplicated = False
    renumeration = []

    for i2 in xrange(0,number_of_nodes2):
        for iRenum in xrange(0, n_duplicated):
            if (i2 == duplicateIndices[iRenum][0]):
                renumeration.append(duplicateIndices[iRenum][1])
                n_renumerated = n_renumerated + 1
                duplicated = True
                break

        if (duplicated == False):
            renumeration.append(i2 + number_of_nodes1 - n_renumerated)
            for component in xrange(0,3):
                vertices_merged[component].append(vertices_two[component][i2])

        duplicated = False

    # merge elements and domain indices
    elements_merged = elements_one.tolist()
    domain_indices_merged = domain_indices_one
    domain_index_addition = max(domain_indices_one) + 1

    number_of_elements2 = len(elements_two[0] + domain_index_addition);

    for i2 in xrange(0,number_of_elements2):
        domain_indices_merged.append(domain_indices_two[i2] + 1)
        for component in xrange(0,3):
            elements_merged[component].append(renumeration[elements_two[component][i2]])

    #create mesh from vertices, elements, domain indices

    from bempp.api.grid import grid_from_element_data

    return grid_from_element_data    (
                                        np.asarray(vertices_merged),
                                        np.asarray(elements_merged),
                                        domain_indices_merged
                                     )

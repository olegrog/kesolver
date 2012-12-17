
# [elm_type, nodes]
facets_of_tetrahedron = [[2, 1, 2, 3], [2, 1, 3, 4], [2, 1, 4, 2], [2, 2, 4, 3]]
facets_of_prism = [[2, 1, 2, 3], [3, 1, 3, 6, 4], [3, 3, 2, 5, 6], [3, 1, 4, 5, 2], [2, 4, 6, 5]]
facets_of_hexahedron = [[3, 1, 2, 3, 4], [3, 1, 5, 6, 2], [3, 2, 6, 7, 3], [3, 3, 7, 8, 4], [3, 1, 4, 8, 5], [3, 5, 8, 7, 6]]

def facets_of_elem(elm):
    dic_of_facets = {4 : facets_of_tetrahedron,
                     5 : facets_of_hexahedron,
                     6 : facets_of_prism}
    facets = dic_of_facets[elm.type]
    for facet in facets:
        yield element(facet[0], [elm.nodes[i-1] for i in facet[1:]])

class element:
    def __init__(self, t, ns, ord_index=0, phys_index=0, part_index=0, neigbors=None):
        self.type = t
        self.nodes = ns
        self.ord_index  = ord_index
        self.phys_index = phys_index
        self.part_index = part_index
        self.neigbors   = neigbors if neigbors is not None else []


    @staticmethod
    def number_of_nodes(t):
        number_of_nodes_list = [-1, 2, 3, 4, 4, 8, 6, 5, 3, 6, 9, 10, 27, 18, 14, 1]
        return number_of_nodes_list[t]
        

    def dimension(self):
        dimensions_list = [-1, 1, 2, 2, 3, 3, 3, 3, 1, 2, 2, 3, 3, 3, 3, 0] 
        return dimensions_list[self.type]


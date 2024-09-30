import numpy as np
from modules.orthonormal import orthonormal_pair, face_vector
import numpy.linalg as la
from numpy.linalg import norm

class layer:
    def __init__(self, monomer_diameter, start_pos, heading, starting_index:int):
        self.__monomer_diameter = monomer_diameter
        self.__start_pos = np.array(start_pos)
        self.__heading = np.array(heading)
        
        self.positions = []
        self.positions.append(self.__start_pos)
        
        self.__start_index = int(starting_index)
        self.indices = []
        self.indices.append(self.__start_index)
        
        self.__generate_layer()
        
    
    def __generate_basis(self):
        h = self.__heading
        h = h / norm(h)
        
        f, g = orthonormal_pair(h)
        
        return h, f, g
    
    def __generate_layer(self):
        h, f, g = self.__generate_basis()
        
        a = self.__monomer_diameter
        
        # +g => 2;  +(f+g) => 3; +f => 4
        p1 = self.__start_pos
        p2 = p1 + a*g
        p3 = p1 + a* (f+g)
        p4 = p1 + a*f
        
        index1 = int(self.__start_index)
        index2 = int(index1 + 1)
        index3 = int(index1 + 2)
        index4 = int(index1 + 3)
        
        # p1 is already appended to self.positions when creating the layer, start from p2
        self.positions.append(p2)
        self.positions.append(p3)
        self.positions.append(p4)
        
        # index1 is already appended to self.indices when creating the layer, start from index2
        self.indices.append(index2)
        self.indices.append(index3)
        self.indices.append(index4)
    
    def make_next_layer(self):
        a = self.__monomer_diameter
        h = self.__heading
        h = h / norm(h)
        start_pos = self.positions[0] + a * h
        start_index = int(self.indices[0] + 4)
        return layer(self.__monomer_diameter, start_pos, self.__heading, start_index)

    def get_basis(self):
        h, f, g = self.__generate_basis()
        return h, f, g
    
    def get_indices(self):
        return self.indices

    @property
    def mono_diameter(self):
        return self.__monomer_diameter
    
    @property
    def atom_types(self):
        return self.__atom_types
        
        
        
        
        
        
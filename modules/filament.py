import numpy as np
from modules.layer import layer
import numpy.linalg as la
from numpy.linalg import norm
from numpy import sqrt, pi, cos, sin, arcsin, arccos


class filament:
    def __init__(self, num_monomers, monomer_diameter, start_pos, heading, linker_distance, linker_list, linker_diameter, radius_of_curvature):
        self.__num_monomers = num_monomers
        self.__monomer_diameter = monomer_diameter

        self.__start_pos = np.array(start_pos)
        self.__heading = np.array(heading)

        self._radius_of_curvature = radius_of_curvature

        self.__layers = []
        self.__monomer_layer_units = []

        self.__basis_sets = []
        self.__angles = []
        self.__bonds = []

        self.__linker_list = linker_list
        self.__linker_distance = linker_distance
        self.__linker_diameter = linker_diameter
        self.__num_linkers = sum(linker_list)
        self.__linkers = []
        self.linker_positions = []

        try:
            assert len(self.__linker_list) == self.__num_monomers
        except AssertionError:
            print("Error: Must specify whether each monomer has a linker or not")
            print("Number of monomers: ", self.__num_monomers)
            print("Number of linkers specified: ", len(self.__linker_list))
            raise

        self.__generate_filament()

    def __generate_filament(self):
        # Monomer circle of diameter a is inscribed in a square of side length a
        a = (self.__monomer_diameter)

        starting_index = 1

        this_layer = layer(a, self.__start_pos, self.__heading, starting_index)
        self.__layers.append(this_layer)
        h, f, g = this_layer.get_basis()
        self.__basis_sets.append(h)
        self.__basis_sets.append(f)
        self.__basis_sets.append(g)

        # Generating layers to build the filament
        num_layers = len(self.__layers)
        for _ in range(1, self.__num_monomers + 1):
            this_layer = this_layer.make_next_layer()
            self.__layers.append(this_layer)
            num_layers += 1

        # Create monomers as association of 2 layers
        monomer_index = 0
        for i in range(num_layers - 1):
            has_linker = self.__linker_list[monomer_index]
            unit = [monomer_index, i, i+1, has_linker]
            self.__monomer_layer_units.append(unit)
            monomer_index += 1
        try:
            assert len(self.__monomer_layer_units) == self.__num_monomers
        except AssertionError:
            print("Error: Number of monomers do not match number of monomer layer units")
            print("Number of monomers: ", self.__num_monomers)
            print("Number of monomer layer units: ",
                  len(self.__monomer_layer_units))
            raise

        # FLAT LINKERS MODEL
        # Create linker positions

        R = self._radius_of_curvature
        d = self.__linker_distance

        a1 = a + ((a**2)) / np.sqrt(4 * R**2 - a**2)
        a2 = a - ((a**2)) / np.sqrt(4 * R**2 - a**2)

        l = a * R / np.sqrt(4 * R**2 - a**2)

        aF = a1
        aL = self.__linker_diameter

        s1 = (aF - aL) / 2
        s2 = (1/2) * np.sqrt((aF - aL)**2 + d**2)

        theta1 = arcsin(a / (2 * l))
        theta2 = pi - theta1

        gamma = sqrt((1/4)*(a - aL)**2 + (1/4) * (a1 - aL)**2 + d**2)

        phi1 = arccos((a1 - aL) / (2 * gamma))
        phi2 = pi - phi1

        phi3 = arccos((a - aL) / (2 * gamma))
        phi4 = pi - phi3

        # Generate layers for linkers
        for unit_i, unit in enumerate(self.__monomer_layer_units):
            monomer_index, layer_i, next_layer_i, has_linker = unit
            if has_linker:
                p2 = self.__layers[layer_i].positions[1]

                linker_heading = g
                linker_starting_pos = p2 + d * g + abs(s1) * (h+f)

                linker_starting_index = (
                    len(self.__layers) * 4) + (sum(self.__linker_list[:monomer_index]) * 4) + 1

                linker = layer(self.__linker_diameter, linker_starting_pos,
                               linker_heading, linker_starting_index)
                self.__linkers.append(linker)

        # Save positions for each point on linker
        for linker_i in range(len(self.__linkers)):
            linker = self.__linkers[linker_i]
            self.linker_positions.append(linker.positions)

        # Create bonds and angles

        linker_index = 0  # index of the linker in the list of linkers
        # This number may be less than the number of monomers, so must be incremented by one iff has_linker is True

        for unit_i, unit in enumerate(self.__monomer_layer_units):
            monomer_index, layer_i, next_layer_i, has_linker = unit

            layer1 = self.__layers[layer_i]
            layer2 = self.__layers[next_layer_i]

            i1 = int(layer1.indices[0])
            i2 = int(layer1.indices[1])
            i3 = int(layer1.indices[2])
            i4 = int(layer1.indices[3])
            i5 = int(layer2.indices[0])
            i6 = int(layer2.indices[1])
            i7 = int(layer2.indices[2])
            i8 = int(layer2.indices[3])

            # Creating bond connections for the cube when linkers are not involved

            # Intra-layer bonds
            bond_type = 1
            self.__bonds.append([bond_type, i1, i4])
            self.__bonds.append([bond_type, i2, i3])
            if monomer_index == self.__num_monomers - 1:
                self.__bonds.append([bond_type, i5, i8])
                self.__bonds.append([bond_type, i6, i7])

            bond_type = 4
            self.__bonds.append([bond_type, i1, i2])
            self.__bonds.append([bond_type, i3, i4])
            if monomer_index == self.__num_monomers - 1:
                self.__bonds.append([bond_type, i5, i6])
                self.__bonds.append([bond_type, i7, i8])

            # Inter-layer bonds
            bond_type = 2
            self.__bonds.append([bond_type, i1, i5])
            self.__bonds.append([bond_type, i4, i8])

            bond_type = 3
            self.__bonds.append([bond_type, i2, i6])
            self.__bonds.append([bond_type, i3, i7])

            # Create angles when linkers are not involved
            angle_type = 1  # theta1

            self.__angles.append([angle_type, i1, i2, i6])
            self.__angles.append([angle_type, i5, i6, i2])
            self.__angles.append([angle_type, i4, i3, i7])
            self.__angles.append([angle_type, i8, i7, i3])

            # _____________________________________________

            angle_type = 2  # theta2

            self.__angles.append([angle_type, i2, i1, i5])
            self.__angles.append([angle_type, i6, i5, i1])
            self.__angles.append([angle_type, i3, i4, i8])
            self.__angles.append([angle_type, i7, i8, i4])

            # _____________________________________________

            angle_type = 3  # pi/2

            self.__angles.append([angle_type, i1, i2, i3])
            self.__angles.append([angle_type, i2, i3, i4])
            self.__angles.append([angle_type, i3, i4, i1])
            self.__angles.append([angle_type, i4, i1, i2])
            self.__angles.append([angle_type, i5, i6, i7])
            self.__angles.append([angle_type, i6, i7, i8])
            self.__angles.append([angle_type, i7, i8, i5])
            self.__angles.append([angle_type, i8, i5, i6])

            # ---------------------------------------------

            self.__angles.append([angle_type, i1, i5, i8])
            self.__angles.append([angle_type, i5, i8, i4])
            self.__angles.append([angle_type, i8, i4, i1])
            self.__angles.append([angle_type, i4, i1, i5])

            # ---------------------------------------------

            self.__angles.append([angle_type, i2, i6, i7])
            self.__angles.append([angle_type, i6, i7, i3])
            self.__angles.append([angle_type, i7, i3, i2])
            self.__angles.append([angle_type, i3, i2, i6])

            # Bonds and angles involving linkers

            # linker_index = 0  # index of the linker in the list of linkers
            # # This number may be less than the number of monomers, so must be incremented by one iff has_linker is True

            if has_linker:
                linker = self.__linkers[linker_index]
                linker_indices = linker.indices
                l1 = int(linker_indices[0])
                l2 = int(linker_indices[1])
                l3 = int(linker_indices[2])
                l4 = int(linker_indices[3])

                linker_index += 1  # since has_linker is True

                # intra-layer bonds
                bond_type = 6
                self.__bonds.append([bond_type, l1, l2])
                self.__bonds.append([bond_type, l2, l3])
                self.__bonds.append([bond_type, l3, l4])
                self.__bonds.append([bond_type, l4, l1])

                # linker-monomer bonds
                bond_type = 5
                self.__bonds.append([bond_type, i2, l1])
                self.__bonds.append([bond_type, i6, l2])
                self.__bonds.append([bond_type, i7, l3])
                self.__bonds.append([bond_type, i3, l4])

                # Create angles when linkers are involved

                angle_type = 3  # pi/2

                self.__angles.append([angle_type, l1, l2, l3])
                self.__angles.append([angle_type, l2, l3, l4])
                self.__angles.append([angle_type, l3, l4, l1])
                self.__angles.append([angle_type, l4, l1, l2])

                # _____________________________________________

                angle_type = 4  # phi1

                self.__angles.append([angle_type, l1, i2, i6])
                self.__angles.append([angle_type, l2, i6, i2])
                self.__angles.append([angle_type, l4, i3, i7])
                self.__angles.append([angle_type, l3, i7, i3])

                # _____________________________________________

                angle_type = 5  # phi2

                self.__angles.append([angle_type, i2, l1, l2])
                self.__angles.append([angle_type, i6, l2, l1])
                self.__angles.append([angle_type, i3, l4, l3])
                self.__angles.append([angle_type, i7, l3, l4])

                # _____________________________________________

                angle_type = 6  # phi3

                self.__angles.append([angle_type, l2, i6, i7])
                self.__angles.append([angle_type, l3, i7, i6])
                self.__angles.append([angle_type, l1, i2, i3])
                self.__angles.append([angle_type, l4, i3, i2])

                # _____________________________________________

                angle_type = 7  # phi4

                self.__angles.append([angle_type, i6, l2, l3])
                self.__angles.append([angle_type, i7, l3, l2])
                self.__angles.append([angle_type, i2, l1, l4])
                self.__angles.append([angle_type, i3, l4, l1])

    #############################################################################
    # Methods for the filament
    #############################################################################

    def get_parameters(self):
        a = (self.__monomer_diameter)
        R = self._radius_of_curvature
        d = self.__linker_distance

        a1 = a + ((a**2)) / np.sqrt(4 * R**2 - a**2)
        a2 = a - ((a**2)) / np.sqrt(4 * R**2 - a**2)

        l = a * R / np.sqrt(4 * R**2 - a**2)

        aF = a1
        aL = self.__linker_diameter

        s1 = (aF - aL) / 2
        s2 = (1/2) * np.sqrt((aF - aL)**2 + d**2)

        theta1 = arcsin(a / (2 * l))
        theta2 = pi - theta1

        gamma = sqrt((1/4)*(a - aL)**2 + (1/4) * (a1 - aL)**2 + d**2)

        phi1 = arccos((a1 - aL) / (2 * gamma))
        phi2 = pi - phi1

        phi3 = arccos((a - aL) / (2 * gamma))
        phi4 = pi - phi3

        return R, d, a, a1, a2, l, s1, s2, aF, aL, theta1, theta2, gamma, phi1, phi2, phi3, phi4

    #############################################################################
    # Properties of the filament
    #############################################################################

    @property
    def basis_sets(self):
        return self.__basis_sets

    @property
    def radius_of_curvature(self):
        return self._radius_of_curvature

    @property
    def linker_distance(self):
        return self.__linker_distance

    @property
    def monomer_layer_units(self):
        return self.__monomer_layer_units

    @property
    def layers(self):
        return self.__layers

    @property
    def monomer_diameter(self):
        return self.__monomer_diameter

    @property
    def bonds(self):
        return self.__bonds

    @property
    def angles(self):
        return self.__angles

    @property
    def linker_list(self):
        return self.__linker_list

    @property
    def num_linkers(self):
        return self.__num_linkers

    @property
    def linkers(self):
        return self.__linkers

    @property
    def linker_diameter(self):
        return self.__linker_diameter

    @property
    def num_monomers(self):
        return self.__num_monomers

    @property
    def num_layers(self):
        return len(self.__layers)

    @property
    def num_bonds(self):
        return len(self.__bonds)

    @property
    def num_angles(self):
        return len(self.__angles)

    @property
    def total_particles(self):
        layers_particles = 4 * len(self.__layers)
        linkers_particles = 4 * len(self.__linkers)
        return layers_particles + linkers_particles

    @property
    def heading(self):
        return self.__heading

    @property
    def start_pos(self):
        return self.__start_pos

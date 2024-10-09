import numpy as np
from modules.angle import filter_angle
from modules.filament import filament


def write_polymer_data(filament_name: filament, box_dimensions: list, mass: list, bond_styles:list, angle_styles:list, data_fname_str: str):
    
    R, d, a, a1, a2, l, s1, s2, aF, aL, theta1, theta2, gamma, phi1, phi2, phi3, phi4 = filament_name.get_parameters()

    num_monomers = filament_name.num_monomers
    num_layers = filament_name.num_layers
    num_total_particles = filament_name.total_particles
    num_linkers = filament_name.num_linkers
    num_bonds = filament_name.num_bonds
    num_angles = filament_name.num_angles
    
    num_atom_types = len(mass)
    num_bond_types = len(bond_styles)
    num_angle_types = len(angle_styles)
    
    xlo, xhi, ylo, yhi, zlo, zhi = box_dimensions

    with open(data_fname_str, 'w') as data_f:
        
        # Header
        data_f.write("LAMMPS data for trapezoidal filament with flat linker\n\n")
        
        # Numbers
        data_f.write("{} atoms\n".format(num_total_particles))
        data_f.write("{} bonds\n".format(num_bonds))
        data_f.write("{} angles\n".format(num_angles))
        
        data_f.write("\n")
        
        # Types
        data_f.write("{} atom types\n".format(num_atom_types))
        data_f.write("{} bond types\n".format(num_bond_types))
        data_f.write("{} angle types\n".format(num_angle_types))
        
        data_f.write("\n")
        
        # Box dimensions
        data_f.write("{:.4f} {:.4f} xlo xhi\n".format(xlo, xhi))
        data_f.write("{:.4f} {:.4f} ylo yhi\n".format(ylo, yhi))
        data_f.write("{:.4f} {:.4f} zlo zhi\n".format(zlo, zhi))
        
        data_f.write("\n")
        
        # Masses
        data_f.write("Masses\n\n")
        
        for m in mass:
            particle_type, particle_mass, particle_name = m
            data_f.write("{:d} {:.4f} # {}\n".format(particle_type, particle_mass, particle_name))
        
        data_f.write("\n")
        
        # Atoms
        data_f.write("Atoms\n\n")
        
        molecule_index = 1
        
        atom_type = 1 # Monomer vertex
        for layer_i, layer in enumerate(filament_name.layers):
            for j in range(4):
                atom_index = layer.indices[j]
                px, py, pz = layer.positions[j][0], layer.positions[j][1], layer.positions[j][2]
                data_f.write("{:d} {:d} {:d} {:.4f} {:.4f} {:.4f}\n".format(atom_index, molecule_index, atom_type, px, py, pz))

        atom_type = 2 # Linker vertex
        for linker_i, linker in enumerate(filament_name.linkers):
            for atom_i in range(len(linker.positions)):
                atom_index = linker.indices[atom_i]
                px, py, pz = linker.positions[atom_i]
                data_f.write("{:d} {:d} {:d} {:.4f} {:.4f} {:.4f}\n".format(atom_index, molecule_index, atom_type, px, py, pz))
        
        # for monomer in filament_name.monomer_layer_units:
        #     monomer_index, layer1_i, layer2_i, has_linker = monomer
        #     l1, l2 = filament_name.layers[layer1_i], filament_name.layers[layer2_i]
            
        #     atom_type = 1
            
        #     for atom_i in range(len(l1.positions)):
        #         atom_index = l1.indices[atom_i]
        #         px, py, pz = l1.positions[atom_i]
        #         data_f.write("{:d} {:d} {:d} {:.4f} {:.4f} {:.4f}\n".format(atom_index, molecule_index, atom_type, px, py, pz))
            
        #     if layer2_i == num_monomers:
        #         atom_type = 1
        #         for atom_i in range(len(l2.positions)):
        #             atom_index = l2.indices[atom_i]
        #             px, py, pz = l2.positions[atom_i]
        #             data_f.write("{:d} {:d} {:d} {:.4f} {:.4f} {:.4f}\n".format(atom_index, molecule_index, atom_type, px, py, pz))
            
        #     if has_linker:
        #         atom_type = 2
        #         linker = filament_name.linkers[monomer_index]
        #         for atom_i in range(len(linker.positions)):
        #             atom_index = linker.indices[atom_i]
        #             px, py, pz = linker.positions[atom_i]
        #             data_f.write("{:d} {:d} {:d} {:.4f} {:.4f} {:.4f}\n".format(atom_index, molecule_index, atom_type, px, py, pz))
        
        data_f.write("\n")
        
        # Bonds
        
        data_f.write("Bonds\n\n")
        
        for bond_i, bond in enumerate(filament_name.bonds):
            bond_type, i1, i2 = bond
            data_f.write("{:d} {:d} {:d} {:d}\n".format(bond_i+1, bond_type, i1, i2))
        
        data_f.write("\n")
        
        # Angles
        
        data_f.write("Angles\n\n")
        
        for angle_i, angle in enumerate(filament_name.angles):
            angle_type, i1, i2, i3 = angle
            data_f.write("{:d} {:d} {:d} {:d} {:d}\n".format(angle_i+1, angle_type, i1, i2, i3))
        
        
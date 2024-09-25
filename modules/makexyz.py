import numpy as np
from modules.filament import filament


def dump_filament(filepath, filament_list):
    with open(filepath, 'w') as f:
        fc = 0
        n_atoms = 0
        for fn in filament_list:
            n_atoms += fn.total_particles

        f.write("{}\n".format(n_atoms))
        f.write("Properties=atom_types:I:1:molecule:I:1:pos:R:3:radius:R:1\n")
        for filamentname in filament_list:
            fc += 1
            print("Writing filament {} to {}".format(fc, filepath))
            
            atom_type = "m"
            for layer in filamentname.layers:
                    for j in range(4):
                        f.write("{} {} {:f} {:f} {:f} {}\n".format(atom_type, fc, layer.positions[j][0], layer.positions[j][1], layer.positions[j][2], filamentname.monomer_diameter/5))
            
            atom_type = "l"
            for linker_i, linker in enumerate(filamentname.linkers):
                for atom_i in range(len(linker.positions)):
                    px, py, pz = linker.positions[atom_i]
                    atom_number = linker.indices[atom_i]
                    f.write("{} {} {:f} {:f} {:f} {}\n".format(atom_type, fc, px, py, pz, filamentname.linker_diameter/5))

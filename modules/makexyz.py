import numpy as np
from modules.filament import filament


def dump_filament(filepath, filament_list):
    with open(filepath, 'w') as f:
        filament_index = 0
        n_atoms = 0
        for fn in filament_list:
            n_atoms += fn.total_particles

        f.write("{}\n".format(n_atoms))
        f.write("Properties=atom_types:I:1:molecule:I:1:pos:R:3:radius:R:1\n")
        for filamentname in filament_list:
            filament_index += 1
            print("Writing filament {} to {}".format(filament_index, filepath))
            
            atom_type = "m"
            for layer_i, layer in enumerate(filamentname.layers):
                    for j in range(4):
                        # atom_type = "m{}".format(j+1 + 4*layer_i)
                        px = layer.positions[j][0]
                        py = layer.positions[j][1]
                        pz = layer.positions[j][2]
                        f.write("{} {} {:f} {:f} {:f} {}\n".format(atom_type, filament_index, px, py, pz, filamentname.monomer_diameter/5))
            
            atom_type = "l"
            for linker_i, linker in enumerate(filamentname.linkers):
                for atom_i in range(len(linker.positions)):
                    # atom_type = "l{}".format(atom_i+1)
                    px, py, pz = linker.positions[atom_i]
                    atom_number = linker.indices[atom_i]
                    f.write("{} {} {:f} {:f} {:f} {}\n".format(atom_type, filament_index, px, py, pz, filamentname.linker_diameter/5))

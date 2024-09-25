import numpy as np
from modules.filament import filament

def dump_filament(filepath, filament_list, draw_frame = False):
    with open(filepath, 'w') as f:
        fc = 0
        n_atoms = 0
        for fn in filament_list:
            n_atoms += fn.num_beads
            if draw_frame:
                n_atoms += 4*len(fn.layers)
        f.write("{}\n".format(n_atoms))
        f.write("Properties=atom_types:I:1:molecule:I:1:pos:R:3:radius:R:1\n")
        for filamentname in filament_list:
            fc+=1
            print("Writing filament {} to {}".format(fc, filepath))
            for i in range(filamentname.num_beads):
                print("Writing bead {} of filament {}".format(i, fc))
                if i == 0:
                    atom_type = "s"
                elif i == filamentname.num_beads - 1:
                    atom_type = "e"
                else:
                    atom_type = "m"
                f.write("{} {} {:f} {:f} {:f} {}\n".format(atom_type, fc, filamentname.beads[i][0], filamentname.beads[i][1], filamentname.beads[i][2], filamentname.monomer_diameter/2))
            if draw_frame:
                for layer in filamentname.layers:
                    for j in range(4):
                        atom_type = "f"
                        f.write("{} {} {:f} {:f} {:f} {}\n".format(atom_type, fc, layer.positions[j][0], layer.positions[j][1], layer.positions[j][2], filamentname.monomer_diameter/20))
            
            
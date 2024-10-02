import numpy as np
from modules.angle import filter_angle
from modules.filament import filament

def save_filament_info(filament_name: filament, info_file: str):
    R, d, a, a1, a2, l, s1, s2, aF, aL, theta1, theta2, gamma, phi1, phi2, phi3, phi4 = filament_name.get_parameters()

    num_monomers = filament_name.num_monomers
    num_layers = filament_name.num_layers
    num_total_particles = filament_name.total_particles
    num_linkers = filament_name.num_linkers
    num_bonds = filament_name.num_bonds
    num_angles = filament_name.num_angles
    
    with open('info.txt', 'w') as info_file:
        info_file.write("#R \t d \t a \t a1 \t a2 \t l \t s1 \t s2 \t aF \t aL \t theta1 \t theta2 \t gamma \t phi1 \t phi2 \t phi3 \t phi4 num_monomers \t num_layers \t num_total_particles \t num_linkers \t num_bonds \t num_angles\n")
        
        info_file.write("{:.4f} \t".format(R))
        info_file.write("{:.4f} \t".format(d))
        info_file.write("{:.4f} \t".format(a))
        info_file.write("{:.4f} \t".format(a1))
        info_file.write("{:.4f} \t".format(a2))
        info_file.write("{:.4f} \t".format(l))
        info_file.write("{:.4f} \t".format(s1))
        info_file.write("{:.4f} \t".format(s2))
        info_file.write("{:.4f} \t".format(aF))
        info_file.write("{:.4f} \t".format(aL))
        info_file.write("{:.4f} \t".format(np.degrees(theta1)))
        info_file.write("{:.4f} \t".format(np.degrees(theta2)))
        info_file.write("{:.4f} \t".format(gamma))
        info_file.write("{:.4f} \t".format(np.degrees(phi1)))
        info_file.write("{:.4f} \t".format(np.degrees(phi2)))
        info_file.write("{:.4f} \t".format(np.degrees(phi3)))
        info_file.write("{:.4f} \t".format(np.degrees(phi4)))
        info_file.write("{} \t".format(num_monomers))
        info_file.write("{} \t".format(num_layers))
        info_file.write("{} \t".format(num_total_particles))
        info_file.write("{} \t".format(num_linkers))
        info_file.write("{} \t".format(num_bonds))
        info_file.write("{} \n".format(num_angles))
    
    
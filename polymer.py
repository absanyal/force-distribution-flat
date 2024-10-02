import numpy as np
from modules.angle import filter_angle
from modules.filament import filament
from modules.save_filament_info import save_filament_info
from modules.polymer_data_writer import write_polymer_data

# various toggles
auto_generate_seed = 1
dump_minimization = 1

# Name of polymer data file
data_fname_str = 'polymer.data'

# Name of LAMMPS input file
input_fname_str = 'input.lammps'

# Name of information file
info_fname_str = 'info.txt'

# Name of linker position file
link_pos_fname_str = 'data/link_pos.txt'

# Name of end positions file
end_pos_fname_str = 'data/e2e_pos.txt'

# Name of CoM MSD file
com_msd_fname_str = 'data/com_msd.txt'

# Name of CoM position file
com_pos_fname_str = 'data/com_pos.txt'

# ---Box dimensions---
xlo, xhi = 0.0, 700
ylo, yhi = 0.0, 1000
zlo, zhi = 0.0, 700

# This is passed to write_polymer_data, do not change this line
box_dimensions = [xlo, xhi, ylo, yhi, zlo, zhi]

# ---Filament parameters---
num_monomers = 1
monomer_diameter = 5
linker_distance = 2.5
linker_diameter = 2
radius_of_curvature = 100
distance_from_axis = 325

# Angle of the filament with the wall
angle = 90
angle = np.radians(angle)  # COnverts to radians, do not change this line

# Calculating the start position and heading of the filament
start_pos = [(xhi - xlo)/2.0 + distance_from_axis,
             (yhi - ylo)/2.0, (zhi - zlo)/2.0]
heading = [0, np.cos(angle), -np.sin(angle)]

# Linker list
linker_list = np.ones(num_monomers)

# Create the filament
f1 = filament(num_monomers, monomer_diameter, start_pos, heading,
              linker_distance, linker_list, linker_diameter, radius_of_curvature)

# Get the parameters of the filament for making bonds and angles
R, d, a, a1, a2, l, s1, s2, aF, aL, theta1, theta2, gamma, phi1, phi2, phi3, phi4 = f1.get_parameters()

theta1 = np.degrees(theta1)
theta2 = np.degrees(theta2)
phi1 = np.degrees(phi1)
phi2 = np.degrees(phi2)
phi3 = np.degrees(phi3)
phi4 = np.degrees(phi4)
pi_by_2 = 90

# Particle types
mass = [
    [1, 1.0, "monomer_vertex"],
    [2, 1.0, "linker_vertex"]
]

# Bond styles: Bond type, Bond style, k, r0
bond_styles = [
    [1, "harmonic", 1500.0, a],
    [2, "harmonic", 1500.0, a2],
    [3, "harmonic", 1500.0, a1],
    [4, "harmonic", 1500.0, 2 * l],
    [5, "harmonic", 1500.0, gamma],
    [6, "harmonic", 1500.0, aL]
]

# Angle styles: Angle type, Angle style, k, theta0
angle_styles = [
    [1, "harmonic", 1500.0, theta1],
    [2, "harmonic", 1500.0, theta2],
    [3, "harmonic", 1500.0, pi_by_2],
    [4, "harmonic", 1500.0, phi1],
    [5, "harmonic", 1500.0, phi2],
    [6, "harmonic", 1500.0, phi3],
    [7, "harmonic", 1500.0, phi4]
]

###################################################################################
######################### WRITE POLYMER DATA ######################################
###################################################################################

# ---Info file---
save_filament_info(f1, info_fname_str)

# ---Data file---
write_polymer_data(f1, box_dimensions, mass, bond_styles, angle_styles, data_fname_str)

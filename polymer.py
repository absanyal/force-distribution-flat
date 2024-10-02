import numpy as np
from modules.angle import filter_angle
from modules.filament import filament
from modules.save_filament_info import save_filament_info
from modules.polymer_data_writer import write_polymer_data
from modules.lammps_input_writer import write_lammps_input

###################################################################################
######################### DATA FILE PARAMETERS ####################################
###################################################################################

# various toggles
dump_minimization = 1
auto_generate_seed = 1
fixed_seed = 123456

data_fname_str = 'polymer.data'  # polymer data file
input_fname_str = 'input.lammps'  # lammps input file
info_fname_str = 'info.txt'  # info file
link_pos_fname_str = 'data/link_pos.txt'  # linker positions file
end_pos_fname_str = 'data/e2e_pos.txt'  # end-to-end positions file
com_msd_fname_str = 'data/com_msd.txt'  # CoM MSD file
com_pos_fname_str = 'data/com_pos.txt'  # CoM positions file

###################################################################################
######################### BOX DIMENSIONS ##########################################
###################################################################################

xlo, xhi = 0.0, 700
ylo, yhi = 0.0, 1000
zlo, zhi = 0.0, 700

x_width = xhi - xlo
y_width = yhi - ylo
z_width = zhi - zlo

# This is passed to write_polymer_data, do not change this line
box_dimensions = [xlo, xhi, ylo, yhi, zlo, zhi]

###################################################################################
######################### FILAMENT PARAMETERS #####################################
###################################################################################

num_monomers = 1
monomer_diameter = 5
linker_distance = 2.5
linker_diameter = 2
radius_of_curvature = 100

# Distance of the filament head from the long axis of the cylinder
distance_from_axis = 325

# Angle of the filament with the wall
angle = 90
angle = np.radians(angle)  # Converts to radians, do not change this line

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

# Pair coefficients for hybrid style
# -1 is a placeholder for the pair style zero
pair_coeff = [
    [1, 1, "zero", []],
    [2, 2, "zero", []],
    [1, 2, "lj/cut", [5.0, 2.5, 2.5 * 2.0**(1/6)]]
]

# Pair cutoffs for hybrid style
pair_cutoffs = [
    ["zero", 30.0],
    ["lj/cut", 8.0]
]

# Atom groups
groups = [
    [1, "chain"],
    [2, "linker"]
]

# Wall-atom interactions, all assumed to be LJ type in format epsilon, sigma, cutoff
wall_interactions = [
    ["wallchain", "chain", "lj93", [5.0, 2.5, 2.5]],
    ["walllinker", "linker", "lj93", [5.0, 2.5, 2.5]]
]

###################################################################################
######################### Simulation parameters ###################################
###################################################################################

# Iteration numbers
steps_min = 1000000
steps_run = 1000000

thermo_min = 1000
thermo_run = 1000

record_interval = 1000

dump_interval_min = 1000
dump_interval_run = 1000

temperture = 310.0
time_step = 0.00001

sim_parameters = [steps_min, steps_run, thermo_min, thermo_run,
                  record_interval, dump_interval_min, dump_interval_run, temperture, time_step]

folders = ['data', 'dump', 'link_pos', 'e2e_pos', 'com_pos']

# ----------------------------------------------------------------------------------

# Brownian parameters
if auto_generate_seed:
    seed = np.random.randint(1, 1000000)
else:
    seed = fixed_seed

gamma_t = 1.0

brownian_parameters = [seed, gamma_t]

###################################################################################
######################### WRITE POLYMER DATA ######################################
###################################################################################

# ---Info file---
save_filament_info(f1, info_fname_str)

# ---Data file---
write_polymer_data(f1, box_dimensions, mass, bond_styles,
                   angle_styles, data_fname_str)

###################################################################################
######################### WRITE LAMMPS INPUT ######################################
###################################################################################

# ---LAMMPS input file---
write_lammps_input(filament_name=f1, box_dimensions=box_dimensions, mass=mass, bond_styles=bond_styles, angle_styles=angle_styles, pair_coeff=pair_coeff, pair_cutoffs=pair_cutoffs, groups=groups, wall_interactions=wall_interactions,
                   sim_parameters=sim_parameters, folders=folders, brownian_parameters=brownian_parameters, input_fname_str=input_fname_str, dump_minimization=dump_minimization, filament_datafile=data_fname_str)

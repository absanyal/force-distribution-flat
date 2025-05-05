import numpy as np
from modules.angle import filter_angle
from modules.filament import filament
from modules.save_info import save_filament_info, save_box_info, save_linker_distribution
from modules.polymer_data_writer import write_polymer_data
from modules.lammps_input_writer_langevin import write_lammps_input_langevin

###################################################################################
######################### DATA FILE PARAMETERS ####################################
###################################################################################

# various toggles
dump_minimization = 1
auto_generate_seed = 1
fixed_seed = 123456

data_fname_str = 'polymer.data'  # polymer data file
input_fname_str = 'input.lammps'  # lammps input file
filament_info_fname_str = 'info/filament_info.txt'  # filament info file
box_info_fname_str = 'info/box_info.txt'  # box info file
# linker distribution file
linker_distribution_fname_str = 'info/linker_distribution.txt'

###################################################################################
######################### BOX AND MEMBRANE DIMENSIONS #############################
###################################################################################

xlo, xhi = 0.0, 1000
ylo, yhi = 0.0, 1000
zlo, zhi = 0.0, 1000

x_width = xhi - xlo
y_width = yhi - ylo
z_width = zhi - zlo

# This is passed to write_polymer_data, do not change this line
box_dimensions = [xlo, xhi, ylo, yhi, zlo, zhi]

create_membrane = 0

###################################################################################
######################### FILAMENT PARAMETERS #####################################
###################################################################################

num_monomers = 80
monomer_diameter = 5
linker_distance = 2.5
linker_diameter = 2
radius_of_curvature = 100

# Distance of the filament head from the long axis of the cylinder
# distance_from_axis = 305
distance_from_axis = 0

# Angle of the filament with the wall
angle = 90
angle = np.radians(angle)  # Converts to radians, do not change this line

# Calculating the start position and heading of the filament
start_pos = [(xhi - xlo)/2.0 - distance_from_axis,
             (yhi - ylo)/2.0, (zhi - zlo)/2.0]
heading = [0, np.cos(angle), -np.sin(angle)]

# Linker list
linker_list = np.ones(num_monomers)
# linker_list = np.zeros(num_monomers)
# linker_list = np.random.choice([0, 1], num_monomers)

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
    [1, 1/8, "monomer_vertex"],
    [2, 1/4, "linker_vertex"]
]

# Bond styles: Bond type, Bond style, k, r0
bond_styles = [
    [1, "harmonic", 100.0, a],
    [2, "harmonic", 10.0, a2],
    [3, "harmonic", 10.0, a1],
    [4, "harmonic", 100.0, 2 * l],
    [5, "harmonic", 100.0, gamma],
    [6, "harmonic", 100.0, aL]
]

# Angle styles: Angle type, Angle style, k, theta0
angle_styles = [
    [1, "harmonic", 20, theta1],
    [2, "harmonic", 20, theta2],
    [3, "harmonic", 20, pi_by_2],
    [4, "harmonic", 20, phi1],
    [5, "harmonic", 20, phi2],
    [6, "harmonic", 20, phi3],
    [7, "harmonic", 20, phi4]
]

# Pair coefficients for hybrid style
pair_coeff = [
    [1, 1, "lj/cut", [5.0, 2.5, 2.5 * 2.0**(1/6)]],
    [2, 2, "zero", []],
    [1, 2, "zero", []]
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

###################################################################################
######################### Simulation parameters ###################################
###################################################################################

# Iteration numbers
steps_min = 5000
steps_run = 1000000

thermo_min = 1000
thermo_run = 10000

record_interval = 100

dump_interval_min = 100
dump_interval_run = 1000

temperture = [1.0, 1.0]  # [min, max]
# temperture = [0.0, 0.0]  # [min, max]
time_step = 0.001

# Minimization parameters: [energy_tolerance, force_tolerance, max_iterations, max_evaluations]
energy_tolerance = 0.0
force_tolerance = 1.0e-5
max_iterations = 1000
max_evaluations = 1000
minimization_parameters = [energy_tolerance,
                           force_tolerance, max_iterations, max_evaluations]

sim_parameters = [steps_min, steps_run, thermo_min, thermo_run,
                  record_interval, dump_interval_min, dump_interval_run, temperture, time_step, minimization_parameters]

folders = ['data', 'dump', 'link_pos', 'com_pos', 'mon_pos']

# ----------------------------------------------------------------------------------

# langevin parameters
if auto_generate_seed:
    seed = np.random.randint(1, 1000000)
else:
    seed = fixed_seed

langevin_damp = 1.0  # Langevin damping coefficient

langevin_parameters = [seed, langevin_damp]

# ----------------------------------------------------------------------------------
# Fixes
# Do not include langevin fix for all atoms, it is automatically generated

# Fix 1: nve/limit integration for the minimization
fix_nve_min = ["fix_min", 0.01]
# fix_nve_min = ["fix_min", 0.00001]


# Fix 2: nve/limit integration for the simulation
# fix_nve_run = ["fix_run", 0.01]
fix_nve_run = []

# Fix 3: wall-atom LJ interactions
fix_wall = [
    ["wallchain", "chain", [10.0, 2.1, 2.1 * (2/5)**(1/6)]],
    ["walllinker", "linker", [10.0, 2.1, 3.0]]
]

# ----------------------------------------------------------------------------------
# SHAKE constraints
shake_fix_name = "fixshake"

shake_bonds = [5]
shake_angles = [4, 5, 6]
shake_types = []

# shake_bonds = []
# shake_angles = []
# shake_types = []

shake_tolerance = 0.0001
shake_iterations = 20
shake_print_every = 0

# DO NOT CHANGE THIS LINE
shake_parameters = [shake_fix_name, shake_tolerance, shake_iterations, shake_print_every, shake_bonds, shake_angles, shake_types]

# ----------------------------------------------------------------------------------

###################################################################################
######################### WRITE INFO FILES ########################################
###################################################################################

# ---Filament info file---
save_filament_info(f1, filament_info_fname_str)
save_box_info(box_dimensions, box_info_fname_str)
save_linker_distribution(f1, linker_distribution_fname_str)

###################################################################################
######################### WRITE POLYMER DATA ######################################
###################################################################################

# ---Data file---
write_polymer_data(f1, box_dimensions, mass, bond_styles,
                   angle_styles, data_fname_str)

###################################################################################
######################### WRITE LAMMPS INPUT ######################################
###################################################################################

# ---LAMMPS input file---
write_lammps_input_langevin(filament_name=f1, box_dimensions=box_dimensions, create_membrane=create_membrane, mass=mass, bond_styles=bond_styles, angle_styles=angle_styles, pair_coeff=pair_coeff, pair_cutoffs=pair_cutoffs, groups=groups, sim_parameters=sim_parameters,
                   folders=folders, langevin_parameters=langevin_parameters, input_fname_str=input_fname_str, dump_minimization=dump_minimization, filament_datafile=data_fname_str, fix_nve_min=fix_nve_min, fix_nve_run=fix_nve_run, fix_wall=fix_wall, shake_parameters=shake_parameters)

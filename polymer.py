import numpy as np
from modules.angle import filter_angle
from modules.filament import filament
from modules.save_info import save_filament_info, save_box_info, save_linker_distribution
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
filament_info_fname_str = 'info/filament_info.txt'  # filament info file
box_info_fname_str = 'info/box_info.txt'  # box info file
# linker distribution file
linker_distribution_fname_str = 'info/linker_distribution.txt'

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

num_monomers = 20
monomer_diameter = 5
linker_distance = 2.5
linker_diameter = 2
radius_of_curvature = 100

# Distance of the filament head from the long axis of the cylinder
# distance_from_axis = 330
distance_from_axis = 0

# Angle of the filament with the wall
angle = 90
angle = np.radians(angle)  # Converts to radians, do not change this line

# Calculating the start position and heading of the filament
start_pos = [(xhi - xlo)/2.0 - distance_from_axis,
             (yhi - ylo)/2.0, (zhi - zlo)/2.0]
heading = [0, np.cos(angle), -np.sin(angle)]

# Linker list
# linker_list = np.ones(num_monomers)
# linker_list = np.zeros(num_monomers)
# linker_list = np.random.choice([0, 1], num_monomers)

linker_list = []
for i in range(num_monomers):
    if i % 1 == 0:
        linker_list.append(1)
    else:
        linker_list.append(0)
        
# linker_list = np.zeros(num_monomers)
# gap = 6
# shift = 3
# midpoint = num_monomers // 2
# linker_list[midpoint - shift] = 1
# linker_list[midpoint - shift + gap] = 1

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
    [2, 0.25, "linker_vertex"]
]

# Bond styles: Bond type, Bond style, k, r0
bond_styles = [
    [1, "harmonic", 15000.0, a],
    [2, "harmonic", 15000.0, a2],
    [3, "harmonic", 15000.0, a1],
    [4, "harmonic", 15000.0, 2 * l],
    [5, "harmonic", 15000.0, gamma],
    [6, "harmonic", 15000.0, aL]
]

# Angle styles: Angle type, Angle style, k, theta0
angle_styles = [
    [1, "harmonic", 2500.0, theta1],
    [2, "harmonic", 2500.0, theta2],
    [3, "harmonic", 2500.0, pi_by_2],
    [4, "harmonic", 2500.0, phi1],
    [5, "harmonic", 2500.0, phi2],
    [6, "harmonic", 2500.0, phi3],
    [7, "harmonic", 2500.0, phi4]
]

# Pair coefficients for hybrid style
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

###################################################################################
######################### Simulation parameters ###################################
###################################################################################

# Iteration numbers
steps_min = 2000
steps_run = 100000

thermo_min = 100
thermo_run = 10000

record_interval = 1000

dump_interval_min = 100
dump_interval_run = 1000

temperture = 0.01
time_step = 0.00001

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

# Brownian parameters
if auto_generate_seed:
    seed = np.random.randint(1, 1000000)
else:
    seed = fixed_seed

gamma_t = 1.0

brownian_parameters = [seed, gamma_t]

# ----------------------------------------------------------------------------------
# Fixes
# Do not include Brownian fix for all atoms, it is automatically generated

# Fix 1: nve/limit integration for the minimization
fix_nve_min = ["fix_min", 0.000001]

# Fix 2: wall-atom LJ interactions
fix_wall = [
    ["wallchain", "chain", [5.0, 2.1, 2.1 * 2.0**(1/6)]],
    ["walllinker", "linker", [800.0, 2.1, 2.1 * 2.0**(1/6)]]
]

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
write_lammps_input(filament_name=f1, box_dimensions=box_dimensions, mass=mass, bond_styles=bond_styles, angle_styles=angle_styles, pair_coeff=pair_coeff, pair_cutoffs=pair_cutoffs, groups=groups, sim_parameters=sim_parameters,
                   folders=folders, brownian_parameters=brownian_parameters, input_fname_str=input_fname_str, dump_minimization=dump_minimization, filament_datafile=data_fname_str, fix_nve_min=fix_nve_min, fix_wall=fix_wall)

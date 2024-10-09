import numpy as np
from modules.angle import filter_angle
from modules.filament import filament

def write_lammps_input(filament_name: filament, box_dimensions: list, mass: list, bond_styles:list, angle_styles:list, pair_coeff:list, pair_cutoffs:list, groups:list, wall_interactions:list, sim_parameters:list, folders:list, brownian_parameters:list, input_fname_str: str, filament_datafile:str, dump_minimization: bool, fix_nve_min: list, fix_wall: list):
    
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
    
    steps_min, steps_run, thermo_min, thermo_run, record_interval, dump_interval_min, dump_interval_run, temperture, timestep, minimization_parameters = sim_parameters
    
    brownian_seed, gamma_t = brownian_parameters
    
    with open(input_fname_str, 'w') as input_f:
        
        input_f.write("# LAMMPS input script for trapezoidal filament with flat linker\n\n")
        
        input_f.write("include in.variables\n")

        input_f.write("\n")
        
        input_f.write("variable steps_min equal {:d}\n".format(steps_min))
        input_f.write("variable steps_run equal {:d}\n".format(steps_run))
        input_f.write("variable thermo_min equal {:d}\n".format(thermo_min))
        input_f.write("variable thermo_run equal {:d}\n".format(thermo_run))
        input_f.write("variable record_interval equal {:d}\n".format(record_interval))
        input_f.write("variable dump_interval_min equal {:d}\n".format(dump_interval_min))
        input_f.write("variable dump_interval_run equal {:d}\n".format(dump_interval_run))
        
        input_f.write("\n")
        
        input_f.write("variable temperature equal {:.4f}\n".format(temperture))
        
        input_f.write("\n")
        
        input_f.write("variable brownian_seed equal {:d}\n".format(brownian_seed))
        input_f.write("variable gamma_t equal {:.4f}\n".format(gamma_t))
        
        input_f.write("\n")
        
        input_f.write("units lj\n")
        input_f.write("dimension 3\n")
        input_f.write("neighbor 1.5 bin\n")
        input_f.write("neigh_modify every 1 delay 0 check yes\n")
        input_f.write("boundary p p p\n")
        
        input_f.write("\n")
        
        input_f.write("atom_style molecular\n")
        
        input_f.write("\n")
        
        input_f.write("read_data {}\n".format(filament_datafile))
        
        input_f.write("\n")
        
        #------------------------------------------------------
        
        x_len = xhi - xlo
        y_len = yhi - ylo
        z_len = zhi - zlo
        
        if x_len > y_len and x_len > z_len:
            cyl_dim = "x"
            cyl_c1, cyl_c2 = (yhi + ylo) / 2, (zhi + zlo) / 2
            cyl_rad = y_len / 2
            cyl_lo, cyl_hi = xlo, xhi
        elif y_len > x_len and y_len > z_len:
            cyl_dim = "y"
            cyl_c1, cyl_c2 = (xhi + xlo) / 2, (zhi + zlo) / 2
            cyl_rad = x_len / 2
            cyl_lo, cyl_hi = ylo, yhi
        elif z_len > x_len and z_len > y_len:
            cyl_dim = "z"
            cyl_c1, cyl_c2 = (xhi + xlo) / 2, (yhi + ylo) / 2
            cyl_rad = x_len / 2
            cyl_lo, cyl_hi = zlo, zhi
        else:
            raise ValueError("Provided box dimensions are invalid.")
        
        input_f.write("region membrane cylinder {} {:.4f} {:.4f} {:.4f} {:.4f} {:.4f}\n".format(cyl_dim, cyl_c1, cyl_c2, cyl_rad, cyl_lo, cyl_hi))
        
        input_f.write("\n")
        
        #------------------------------------------------------
        
        # Groups
        for group in groups:
            group_type, group_name = group
            input_f.write("group {} type {}\n".format(group_name, group_type))
        
        input_f.write("\n")
        
        #------------------------------------------------------
        
        # pair styles and coefficients
        input_f.write("pair_style hybrid ")
        for cutoff_val in pair_cutoffs:
            style, cutoff = cutoff_val
            input_f.write("{} {:.4f} ".format(style, cutoff))
        input_f.write("\n")
        
        for coeff_val in pair_coeff:
            i1, i2, pair_style, coeffs = coeff_val
            if pair_style == "zero":
                input_f.write("pair_coeff {} {} {}\n".format(i1, i2, pair_style))
            elif pair_style == "lj/cut":
                epsilon, sigma, Rc = coeffs
                input_f.write("pair_coeff {} {} {} {:.4f} {:.4f} {:.4f}\n".format(i1, i2, pair_style, epsilon, sigma, Rc))
            else:
                raise ValueError("Invalid pair style provided.")
        
        input_f.write("\n")
        
        #------------------------------------------------------
        
        # bond styles
        input_f.write("bond_style harmonic\n")
        for bond_style in bond_styles:
            bond_ID, bond_style, k, r0 = bond_style
            if bond_style == "harmonic":
                input_f.write("bond_coeff {} {:.4f} {:.4f}\n".format(bond_ID, k, r0))
            else:
                raise ValueError("Invalid bond style provided.")
        
        input_f.write("\n")
        
        #------------------------------------------------------
        
        # angle styles
        
        input_f.write("angle_style harmonic\n")
        for angle_style in angle_styles:
            angle_ID, angle_style, k, theta0 = angle_style
            if angle_style == "harmonic":
                input_f.write("angle_coeff {} {:.4f} {:.4f}\n".format(angle_ID, k, theta0))
            else:
                raise ValueError("Invalid angle style provided.")
        
        input_f.write("\n")
        
        #------------------------------------------------------
        
        input_f.write("timestep {:.10f}\n".format(timestep))
        input_f.write("\n")
        
        #------------------------------------------------------
        
        input_f.write("shell mkdir ")
        for folder in folders:
            input_f.write("./{}/ ".format(folder))
        input_f.write("\n")
        
        input_f.write("\n")
        
        #------------------------------------------------------
        
        if dump_minimization:   
            input_f.write("dump minimization all atom ${dump_interval_min} dump/dump.min.${xx}.lammpstrj\n")
            input_f.write("\n")
        
        #------------------------------------------------------
        
        etol, ftol, maxiter, maxeval = minimization_parameters
        input_f.write("minimize {:.6f} {:.6f} {:d} {:d}\n".format(etol, ftol, maxiter, maxeval))
        
        input_f.write("\n")
        
        #------------------------------------------------------
        
        input_f.write("fix {} all nve/limit {:.6f}\n".format(fix_nve_min[0], fix_nve_min[1]))
        
        input_f.write("\n")
        
        for fix_wall_i in fix_wall:
            fix_name, fix_type, fix_params = fix_wall_i
            epsilon, sigma, cutoff = fix_params
            input_f.write("fix {} {} wall/region membrane lj93 {:.4f} {:.4f} {:.4f}\n".format(fix_name, fix_type, epsilon, sigma, cutoff))
        
        input_f.write("\n")
        
        #------------------------------------------------------
        
        input_f.write("thermo_style custom step time temp etotal\n")
        input_f.write("thermo ${thermo_min}\n")
        
        input_f.write("\n")
        
        input_f.write("run ${steps_min}\n")
        
        input_f.write("\n")
        
        input_f.write("unfix {}\n".format(fix_nve_min[0]))
        if dump_minimization:
            input_f.write("undump minimization\n")
        
        input_f.write("\n")
        
        
        #------------------------------------------------------
        
        input_f.write("reset_timestep 0\n")
        input_f.write("\n")
        
        input_f.write("variable tsteps equal time\n")
        input_f.write("\n")
        
        input_f.write("dump dumpall all atom ${dump_interval_run} dump/dump.${xx}.lammpstrj\n")
        input_f.write("\n")
        
        #------------------------------------------------------
        

        #------------------------------------------------------
        
        input_f.write("fix brnfix all brownian ${temperature} ${brownian_seed} gamma_t ${gamma_t}\n")
        
        input_f.write("\n")
        
        #------------------------------------------------------
        
        input_f.write("thermo ${thermo_run}\n")
        input_f.write("run ${steps_run}\n")
        
        input_f.write("\n")
        
        # -----------------------------------------------------
        
        input_f.write("write_data dump/finalstate_hot.lammpstrj\n")
        
        
            
        
        
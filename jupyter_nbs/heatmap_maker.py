import numpy as np
import matplotlib.pyplot as plt

def w_n(n, N, lp, Ebind0, R0, R_cell, r_mono, T=1, kB=1):
    """
    Calculate the weight of a given state with n monomers attached.
    
    Parameters
    ----------
    
    n : int 
        Number of monomers attached to the cell
        
    N : int
        Total number of monomers
    
    lp : float
        Persistence length of the polymer
    
    Ebind0 : float
        Binding energy of a single linker
    
    R0 : float
        Radius of curvature of the filament
    
    R_cell : float
        Radius of the cell
        
    r_mono : float
        Radius of a monomer
        
    T : float
        Temperature, default is 1
        
    kB : float
        Boltzmann constant, default is 1
    
    Returns
    -------
    
    w_n : float
        Weight of the state with n monomers attached
    """
    beta = 1.0 / (kB * T)
    dEbend = r_mono * (n - 1) * (kB * T * lp) * ( (1/R_cell) - (1/R0)  )**2
    dEbind = n * Ebind0
    w_n = (N - n + 1) * np.exp(-beta * (dEbend + dEbind))
    return w_n

def w_n_list(N, lp, Ebind0, R0, R_cell, r_mono, T=1, kB=1):
    """
    Calculate the weight of all states with 1 to N monomers attached.
    
    Parameters
    ----------
    
    N : int
        Total number of monomers
    
    lp : float
        Persistence length of the polymer
    
    Ebind0 : float
        Binding energy of a single linker
    
    R0 : float
        Radius of curvature of the filament
    
    R_cell : float
        Radius of the cell
        
    r_mono : float
        Radius of a monomer
        
    T : float
        Temperature, default is 1
        
    kB : float
        Boltzmann constant, default is 1
    
    Returns
    -------
    
    w_n_list : np.ndarray
        Array of weights of states with 1 to N monomers attached
    """
    n_list = np.arange(1, N+1)
    w_list = np.zeros_like(n_list, dtype=float)
    for n_i, n in enumerate(n_list):
        w_list[n_i] = w_n(n, N, lp, Ebind0, R0, R_cell, r_mono, T, kB)
    return w_list

def Z(N, lp, Ebind0, R0, R_cell, r_mono, T=1, kB=1):
    """
    Calculate the partition function of the system.
    
    Parameters
    ----------
    
    N : int
        Total number of monomers
    
    lp : float
        Persistence length of the polymer
    
    Ebind0 : float
        Binding energy of a single linker
    
    R0 : float
        Radius of curvature of the filament
    
    R_cell : float
        Radius of the cell
        
    r_mono : float
        Radius of a monomer
        
    T : float
        Temperature, default is 1
        
    kB : float
        Boltzmann constant, default is 1
    
    Returns
    -------
    
    Z : float
        Partition function of the system
    """
    w_list = w_n_list(N, lp, Ebind0, R0, R_cell, r_mono, T, kB)
    Z = np.sum(w_list)
    return Z

def p_n(n, N, lp, Ebind0, R0, R_cell, r_mono, T=1, kB=1):
    """
    Calculate the probability of a given state with n monomers attached.
    
    Parameters
    ----------
    
    n : int
        Number of monomers attached to the cell
    
    N : int
        Total number of monomers
    
    lp : float
        Persistence length of the polymer
    
    Ebind0 : float
        Binding energy of a single linker
    
    R0 : float
        Radius of curvature of the filament
    
    R_cell : float
        Radius of the cell
        
    r_mono : float
        Radius of a monomer
        
    T : float
        Temperature, default is 1
        
    kB : float
        Boltzmann constant, default is 1
    
    Returns
    -------
    
    p : float
        Probability of the state with n monomers attached
    """
    w = w_n(n, N, lp, Ebind0, R0, R_cell, r_mono, T, kB)
    Z_val = Z(N, lp, Ebind0, R0, R_cell, r_mono, T, kB)
    p = w / Z_val
    return p

def p_n_list(N, lp, Ebind0, R0, R_cell, r_mono, T=1, kB=1):
    """
    Calculate the probability of a given state with n monomers attached.
    
    Parameters
    ----------
    
    N : int
        Total number of monomers
    
    lp : float
        Persistence length of the polymer
    
    Ebind0 : float
        Binding energy of a single linker
    
    R0 : float
        Radius of curvature of the filament
    
    R_cell : float
        Radius of the cell
        
    r_mono : float
        Radius of a monomer
        
    T : float
        Temperature, default is 1
        
    kB : float
        Boltzmann constant, default is 1
    
    Returns
    -------
    
    p_list : np.ndarray
        Array of probabilities of states with 1 to N monomers attached
    """
    n_list = np.arange(1, N+1)
    p_list = np.zeros_like(n_list, dtype=float)
    for n_i, n in enumerate(n_list):
        p_list[n_i] = p_n(n, N, lp, Ebind0, R0, R_cell, r_mono, T, kB)
    return p_list

def average_n(N, lp, Ebind0, R0, R_cell, r_mono, T=1, kB=1):
    """
    Calculate the probability of a given state with n monomers attached.
    
    Parameters
    ----------
    
    N : int
        Total number of monomers
    
    lp : float
        Persistence length of the polymer
    
    Ebind0 : float
        Binding energy of a single linker
    
    R0 : float
        Radius of curvature of the filament
    
    R_cell : float
        Radius of the cell
        
    r_mono : float
        Length of a monomer
        
    T : float
        Temperature, default is 1
        
    kB : float
        Boltzmann constant, default is 1
    
    Returns
    -------
    
    average_n : float
        Average number of monomers attached to the cell
    """
    n_list = np.arange(1, N+1)
    p_list = p_n_list(N, lp, Ebind0, R0, R_cell, r_mono, T, kB)
    average_n = np.sum(n_list * p_list)
    return average_n

def R0_const_N_const_lp_vs_Eb0(lp_min, lp_max, num_lp, Eb0_min, Eb0_max, num_Eb0, N, R0, R_cell, r_mono, mode = 'number', T=1, kB=1):
    """
    Generates a heatmap of average number of monomers attached to the cell as a function of persistence length and binding energy.
    
    Parameters
    ----------
    
    lp_min, lp_max : float
        Minimum and maximum values of the persistence length
    
    num_lp : int
        Number of points to generate between lp_min and lp_max
    
    Eb0_min, Eb0_max : float
        Minimum and maximum values of the binding energy
        
    num_Eb0 : int
        Number of points to generate between Eb0_min and Eb0_max
        
    N : int
        Total number of monomers
        
    R0 : float
        Radius of curvature of the filament
    
    R_cell : float
        Radius of the cell
    
    r_mono : float
        Radius of a monomer
    
    mode : str
        Mode of the plot.
        
        `number:` Number of monomers attached to the cell.
        
        `fraction:` Fraction of monomers attached to the cell.
        
        `status:` Status of the monomer, 1 if at least two monomers are attached, 0 otherwise.
        
    T : float
        Temperature, default is 1
        
    kB : float
        Boltzmann constant, default is 1
        
    Returns
    -------
    
    average_n_matrix : np.ndarray
        Matrix of average number/fraction/status of monomers attached to the cell
    """
    
    lp_list = np.linspace(lp_min, lp_max, num_lp)
    Eb0_list = np.linspace(Eb0_min, Eb0_max, num_Eb0)
    
    average_n_matrix = np.zeros((num_lp, num_Eb0), dtype=float)
    
    for lp_i, lp in enumerate(lp_list):
        for Eb0_i, Eb0 in enumerate(Eb0_list):
            avg_n = average_n(N, lp, Eb0, R0, R_cell, r_mono, T, kB)
            if mode == 'fraction':
                avg_n = avg_n / N
            elif mode == 'status':
                if avg_n >= 2:
                    avg_n = 1
                else:
                    avg_n = 0
            average_n_matrix[lp_i, Eb0_i] = avg_n
    
    return average_n_matrix

def R0_const_lp_const_N_vs_Eb0(N_min, N_max, num_N, Eb0_min, Eb0_max, num_Eb0, lp, R0, R_cell, r_mono, mode = 'number', T=1, kB=1):
    """
    Generates a heatmap of average number of monomers attached to the cell as a function of number of monomers and binding energy.
    
    Parameters
    ----------
    
    N_min, N_max : int
        Minimum and maximum values of the total number of monomers
    
    num_N : int
        Number of points to generate between N_min and N_max
    
    Eb0_min, Eb0_max : float
        Minimum and maximum values of the binding energy
        
    num_Eb0 : int
        Number of points to generate between Eb0_min and Eb0_max
        
    lp : float
        Persistence length of the polymer
        
    R_cell : float
        Radius of the cell
    
    r_mono : float
        Radius of a monomer
    
    mode : str
        Mode of the plot.
        
        `number:` Number of monomers attached to the cell.
        
        `fraction:` Fraction of monomers attached to the cell.
        
        `status:` Status of the monomer, 1 if at least two monomers are attached, 0 otherwise.
        
    T : float
        Temperature, default is 1
        
    kB : float
        Boltzmann constant, default is 1
        
    Returns
    -------
    
    average_n_matrix : np.ndarray
        Matrix of average number/fraction/status of monomers attached to the cell
    """
    
    N_list = np.linspace(N_min, N_max, num_N)
    Eb0_list = np.linspace(Eb0_min, Eb0_max, num_Eb0)
    
    average_n_matrix = np.zeros((num_N, num_Eb0), dtype=float)
    
    for N_i, N in enumerate(N_list):
        for Eb0_i, Eb0 in enumerate(Eb0_list):
            avg_n = average_n(N, lp, Eb0, R0, R_cell, r_mono, T, kB)
            if mode == 'fraction':
                avg_n = avg_n / N
            elif mode == 'status':
                if avg_n >= 2:
                    avg_n = 1
                else:
                    avg_n = 0
            average_n_matrix[N_i, Eb0_i] = avg_n
    
    return average_n_matrix


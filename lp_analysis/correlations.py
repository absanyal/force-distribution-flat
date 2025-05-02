import numpy as np
import matplotlib.pyplot as plt

raw_data = np.loadtxt('../mon_pos/mon_pos.1.txt')

num_time_steps = raw_data.shape[0]
print("Number of time steps: {}".format(num_time_steps))

fraction_of_recording_to_skip = 0.1
t_start_index = int(num_time_steps * fraction_of_recording_to_skip)
print("Starting averaging after {} time steps".format(t_start_index))

t_list = np.arange(1, num_time_steps + 1)
# print(t_list)

t_list = t_list[t_start_index:]
num_time_steps = len(t_list)
print("Number of time steps after skipping: {}".format(num_time_steps))

num_atoms = int((raw_data.shape[1] - 1) / 3)
print("Number of atoms: {}".format(num_atoms))

num_bond_vectors = num_atoms - 1
print("Number of bond vectors: {}".format(num_bond_vectors))

processed_data = np.zeros((num_time_steps, num_bond_vectors, 3))

for t_i in range(num_time_steps):
    row = raw_data[t_i, 1:]
    for a_i in range(num_atoms - 1):
        atom_1_x = row[a_i * 3]
        atom_1_y = row[a_i * 3 + 1]
        atom_1_z = row[a_i * 3 + 2]
        atom_2_x = row[(a_i + 1) * 3]
        atom_2_y = row[(a_i + 1) * 3 + 1]
        atom_2_z = row[(a_i + 1) * 3 + 2]
        bond_vector_x = atom_2_x - atom_1_x
        bond_vector_y = atom_2_y - atom_1_y
        bond_vector_z = atom_2_z - atom_1_z
        bond_vector = [bond_vector_x, bond_vector_y, bond_vector_z]
        bond_vector = np.array(bond_vector)
        
        processed_data[t_i, a_i] = bond_vector
        

s_list = np.arange(0, num_bond_vectors)
# print("s_list: {}".format(s_list))

time_averaged_correlations = np.zeros(len(s_list))

max_t = np.max(t_list)
t_skip = np.round(0.1 * max_t, 0)

for t_i in range(num_time_steps):
    percentage = (t_i + 1) / num_time_steps
    if t_i % t_skip == 0 or t_i == num_time_steps - 1:
        print("Processing step {}/{} | {:.0f}%".format(t_i + 1, num_time_steps, percentage * 100))
    for b_i in range(num_bond_vectors):
        delta_s_range = np.arange(num_bond_vectors - b_i)
        # print("For bond vector {}, correlation will be calculated for bonds:".format(b_i))
        # for delta_s in delta_s_range:
        #     print("Bond vector {}".format(b_i + delta_s))
        for delta_s in delta_s_range:
            v1 = processed_data[t_i, b_i]
            v2 = processed_data[t_i, b_i + delta_s]
            
            v_1_norm = np.linalg.norm(v1)
            v_2_norm = np.linalg.norm(v2)

            corr = np.dot(v1, v2) / (v_1_norm * v_2_norm)
            
            time_averaged_correlations[delta_s] += corr / (num_bond_vectors - delta_s)

# Normalize the time-averaged correlations
time_averaged_correlations /= num_time_steps

with open("lp_data.dat", 'w') as f:
    f.write("# s \t C(s)\n")
    for s_i in range(len(s_list)):
        f.write("{} \t {}\n".format(s_list[s_i], time_averaged_correlations[s_i]))
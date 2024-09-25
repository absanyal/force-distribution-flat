from modules.filament import filament
import numpy as np

num_monomers = 20
monomer_diameter = 5
start_pos = [0,0,0]
heading = [1,0,0]
linker_distance = 2.5
linker_diameter = 5
radius_of_curvature = 100

linker_list = np.random.choice([0,1], num_monomers)
# linker_list = np.ones(num_monomers)
# linker_list = np.zeros(num_monomers)
# linker_list = np.array([0, 1, 0, 1])

f1 = filament(num_monomers, monomer_diameter, start_pos, heading, linker_distance, linker_list, linker_diameter, radius_of_curvature)

print("="*20)

print("Number of monomers: {}".format(f1.num_monomers))
print("Number of layers: {}".format(f1.num_layers))
print("Number of linkers: {}".format(f1.num_linkers))
print("Number of bonds: {}".format(f1.num_bonds))
print("Number of angles: {}".format(f1.num_angles))

print("="*20)

for monomer in f1.monomer_layer_units:
    monomer_index, layer1, layer2, has_linker = monomer
    print("Monomer {} -> Layer {} & Layer {} ; Linker: {}".format(monomer_index, layer1, layer2, int(has_linker)))

print("="*20)

for monomer in f1.monomer_layer_units:
    monomer_index, layer1_i, layer2_i, has_linker = monomer
    l1 = f1.layers[layer1_i]
    l2 = f1.layers[layer2_i]
    print("Indices of layer {} are: {}".format(layer1_i, l1.indices))
    print("Indices of layer {} are: {}".format(layer2_i, l2.indices))

print("="*20)

print("Indices of the linkers are:")
for linker_i, linker in enumerate(f1.linkers):
    print("Linker {} has indices: {}".format(linker_i, linker.indices))

print("="*20)

for bond in f1.bonds:
    bond_type, i1, i2 = bond
    if 1:
        print("Bond type {} between {} and {}".format(bond_type, i1, i2))

print("="*20)

for angle in f1.angles:
    angle_type, i1, i2, i3 = angle
    if 1:
        print("Angle type {} between {}, {} and {}".format(angle_type, i1, i2, i3))

print("="*20)


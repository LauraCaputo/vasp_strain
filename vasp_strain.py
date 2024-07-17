import numpy as np

def read_poscar(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    comment = lines[0].strip()
    scaling_factor = float(lines[1].strip())
    lattice_vectors = np.array([list(map(float, line.split())) for line in lines[2:5]])
    elements = lines[5].split()
    num_atoms = list(map(int, lines[6].split()))
    total_atoms = sum(num_atoms)
    selective_dynamics = False
    start_index = 7
    if lines[7].strip().lower()[0] == 's':
        selective_dynamics = True
        start_index += 1
    if lines[start_index].strip().lower()[0] == 'c':
        cartesian = True
    else:
        cartesian = False
    start_index += 1
    atomic_positions = [line.split() for line in lines[start_index:start_index + total_atoms]]
    
    return (comment, scaling_factor, lattice_vectors, elements, num_atoms, selective_dynamics, cartesian, atomic_positions)

def write_poscar(filename, poscar_data):
    (comment, scaling_factor, lattice_vectors, elements, num_atoms, selective_dynamics, cartesian, atomic_positions) = poscar_data
    
    with open(filename, 'w') as f:
        f.write(f"{comment}\n")
        f.write(f"{scaling_factor:.16f}\n")
        for vec in lattice_vectors:
            f.write("  ".join(f"{x:.16f}" for x in vec) + "\n")
        f.write("  ".join(elements) + "\n")
        f.write("  ".join(map(str, num_atoms)) + "\n")
        if selective_dynamics:
            f.write("Selective dynamics\n")
        if cartesian:
            f.write("Cartesian\n")
        else:
            f.write("Direct\n")
        for pos in atomic_positions:
            f.write("  ".join(str(x) for x in pos) + "\n")

def apply_strain(lattice_vectors, strain_tensor):
    return np.dot(lattice_vectors, strain_tensor)

def adjust_positions(atomic_positions, old_lattice, new_lattice, cartesian):
    old_lattice_inv = np.linalg.inv(old_lattice)
    
    adjusted_positions = []
    for pos in atomic_positions:
        position = np.array(list(map(float, pos[:3])))
        if cartesian:
            fractional_position = np.dot(old_lattice_inv, position)
        else:
            fractional_position = position
        
        new_position = np.dot(new_lattice, fractional_position)
        if not cartesian:
            new_position = fractional_position
        
        adjusted_pos = list(new_position) + pos[3:]
        adjusted_positions.append(adjusted_pos)
    
    return adjusted_positions

def main(poscar_in, poscar_out, strain_tensor):
    poscar_data = read_poscar(poscar_in)
    comment, scaling_factor, lattice_vectors, elements, num_atoms, selective_dynamics, cartesian, atomic_positions = poscar_data
    
    strained_lattice_vectors = apply_strain(lattice_vectors, strain_tensor)
    adjusted_atomic_positions = adjust_positions(atomic_positions, lattice_vectors, strained_lattice_vectors, cartesian)
    
    strained_poscar_data = (comment, scaling_factor, strained_lattice_vectors, elements, num_atoms, selective_dynamics, cartesian, adjusted_atomic_positions)
    write_poscar(poscar_out, strained_poscar_data)
    print(f"Strained POSCAR written to {poscar_out}")

if __name__ == "__main__":
    poscar_in = input("Enter the name of the input POSCAR file: ")
    poscar_out = input("Enter the name of the output POSCAR file: ")
    

    strain_a = float(input("Enter strain percentage for a-axis:")) / 100.0
    strain_b = float(input("Enter strain percentage for b-axis:")) / 100.0
    strain_c = float(input("Enter strain percentage for c-axis:")) / 100.0

    shear_ab = float(input("Enter shear strain percentage for ab plane:")) / 100.0
    shear_bc = float(input("Enter shear strain percentage for bc plane:")) / 100.0
    shear_ca = float(input("Enter shear strain percentage for ca plane:")) / 100.0
    

    strain_tensor = np.array([[1 + strain_a, shear_ab, shear_ca], 
                              [shear_ab, 1 + strain_b, shear_bc], 
                              [shear_ca, shear_bc, 1 + strain_c]])
    
    main(poscar_in, poscar_out, strain_tensor)

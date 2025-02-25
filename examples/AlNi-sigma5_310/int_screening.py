import subprocess
import numpy as np
import pandas as pd

## define varibales

original_data = '5_1.data'
GB_atom_type = 1  ## int
min_volume = 10.8  ## folat
min_distance = 2.0  ## folat


# Step 10: Run LAMMPS script
def run_lammps(lammps_executable="lmp", input_script="in.min"):
    try:
        result = subprocess.run([lammps_executable, "-in", input_script],
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print("LAMMPS output:\n", result.stdout)
        if result.stderr:
            print("LAMMPS errors:\n", result.stderr)
    except FileNotFoundError:
        print(f"Error: The executable '{lammps_executable}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Step 11: Extract GB atoms
def extract_gb_atoms(input_file='9_combined_structure.data', output_file='11_GB_atoms.txt'):
    positions = []
    with open(input_file, 'r') as file:
        for _ in range(17):
            next(file)  # Skip the first 17 lines
        for line in file:
            parts = line.split()
            if len(parts) >= 5 and int(parts[1]) == GB_atom_type:  # Check for type 1 atoms
                x, y, z = parts[2:5]
                positions.append(f"{x} {y} {z}")
    with open(output_file, 'w') as output:
        output.write("\n".join(positions))
    print(f'GB atoms extracted to {output_file}')

# Step 12: Extract interstitial sites with volume > 10
def extract_interstitial_sites(input_file='0.dump', output_file='12_interstitial_atoms.txt'):
    positions = []
    with open(input_file, 'r') as file:
        for _ in range(9):
            next(file)  # Skip the first 9 lines
        for line in file:
            parts = line.split()
            if len(parts) >= 6 and float(parts[5]) > min_volume:  # Check for volume > 10
                x, y, z = parts[2:5]
                positions.append(f"{x} {y} {z}")
    with open(output_file, 'w') as output:
        output.write("\n".join(positions))
    print(f'Interstitial atoms with volume > 10 saved to {output_file}')

# Step 13: Calculate minimum distances between interstitial sites and GB atoms
def calculate_min_distances(int_sites_file='12_interstitial_atoms.txt', gb_atoms_file='11_GB_atoms.txt', output_file='13_min_distances_chunked.txt'):
    A_coords = np.loadtxt(int_sites_file)
    B_coords = np.loadtxt(gb_atoms_file)
    distances = []

    chunk_size = 1000  # Process in chunks to save memory
    for start in range(0, len(A_coords), chunk_size):
        A_chunk = A_coords[start:start + chunk_size]
        diffs = A_chunk[:, np.newaxis, :] - B_coords[np.newaxis, :, :]
        squared_diffs = np.sum(diffs ** 2, axis=-1)
        min_distances = np.sqrt(squared_diffs).min(axis=1)
        distances.extend(min_distances)

    combined_data = np.hstack((A_coords, np.array(distances).reshape(-1, 1)))
    np.savetxt(output_file, combined_data, fmt="%.10f", header="x y z min_distance")
    print(f"Minimum distances calculated and saved to {output_file}")

# Step 14: Remove interstitial sites with min distance < 2
def filter_interstitial_sites(input_file='13_min_distances_chunked.txt', output_file='14_final_int.txt'):
    filtered_positions = []
    data = np.loadtxt(input_file)
    for row in data:
        if row[3] > min_distance:  # Min distance check
            filtered_positions.append(f"{row[0]} {row[1]} {row[2]}")
    with open(output_file, 'w') as output:
        output.write("\n".join(filtered_positions))
    print(f"Filtered interstitial atoms saved to {output_file}")

# Step 15: Add ID and type to interstitial sites
def add_id_type(input_file='14_final_int.txt', reference_file='Al.txt', output_file='15_modified_int_sites.txt', atom_type=3):
    existing_data = np.loadtxt(reference_file)
    starting_id = existing_data.shape[0] + 1  # Start ID after existing atoms
    data = np.loadtxt(input_file)
    num_rows = data.shape[0]

    row_numbers = np.arange(starting_id, starting_id + num_rows).reshape(-1, 1)
    type_column = np.full((num_rows, 1), atom_type)
    modified_data = np.hstack((row_numbers, type_column, data))
    np.savetxt(output_file, modified_data, fmt='%d %d %.10f %.10f %.10f')
    print(f"IDs and types added to interstitial sites, saved to {output_file}")

# Step 16: Modify and sort original data by type
def modify_and_sort(input_file= original_data, output_file='16_modified_Al.data'):
    with open(input_file, 'r') as file:
        header_lines = [next(file) for _ in range(11)]
    df = pd.read_csv(input_file, delim_whitespace=True, header=None, skiprows=11).iloc[:, :8]
    df = df.sort_values(by=df.columns[1]).reset_index(drop=True)
    df.iloc[:, 0] = range(1, len(df) + 1)  # Reindex atom IDs
    with open(output_file, 'w') as file:
        file.writelines(header_lines)
        df.to_csv(file, sep='\t', header=False, index=False)
    print(f"Modified and sorted data saved to {output_file}")

# Step 17: Output the final NC structure
def output_final_structure(al_data='16_modified_Al.data', int_sites='15_modified_int_sites.txt', output_file='17_final_atomic_structure.data'):
    with open(al_data, 'r') as f:
        original_data = f.readlines()
    with open(int_sites, 'r') as f:
        int_data = f.readlines()

    total_atoms = int(int_data[-1].split()[0])
    combined_data = original_data + int_data
    combined_data[2] = f"{total_atoms} atoms\n"
    combined_data[3] = "3 atom types\n"
    masses_section = """Masses

1 26.9815  # Al
2 26.9815  # Al
3 58.74  # Ni\n\n"""
    combined_data.insert(9, masses_section)
    with open(output_file, 'w') as f:
        f.writelines(combined_data)
    print(f"Combined data file created as '{output_file}'")

# Run all steps in order
run_lammps()
extract_gb_atoms()
extract_interstitial_sites()
calculate_min_distances()
filter_interstitial_sites()
add_id_type()
modify_and_sort()
output_final_structure()

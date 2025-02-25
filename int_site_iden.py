import numpy as np
import pandas as pd
from scipy.spatial import Voronoi
from sklearn.cluster import DBSCAN

# Step 1: Load atom positions and determine Voronoi vertices
filename = 'Al.txt'
data = pd.read_csv(filename, delim_whitespace=True, header=None, names=['atom_id', 'x', 'y', 'z'])

# Calculate the min and max for each coordinate axis
xmin_data, xmax_data = data['x'].min(), data['x'].max()
ymin_data, ymax_data = data['y'].min(), data['y'].max()
zmin_data, zmax_data = data['z'].min(), data['z'].max()

# Load atomic coordinates (skipping the atom ID column)
coordinates = np.loadtxt('Al.txt', usecols=(1, 2, 3))
vor = Voronoi(coordinates)

# Dictionary to store vertices for each atom
atom_vertices = {}

# Check if a point is within bounding box limits
def within_bounds(vertex, xmin, xmax, ymin, ymax, zmin, zmax):
    return xmin <= vertex[0] <= xmax and ymin <= vertex[1] <= ymax and zmin <= vertex[2] <= zmax

# Collect vertices for each atom
for i, point in enumerate(vor.points):
    region_index = vor.point_region[i]
    vertices = vor.regions[region_index]
    finite_vertices = [
        vor.vertices[v].tolist() for v in vertices
        if v != -1 and within_bounds(vor.vertices[v], xmin_data, xmax_data, ymin_data, ymax_data, zmin_data, zmax_data)
    ]
    atom_vertices[i] = finite_vertices

# Save vertices to file
with open('1_Al_vertice.txt', 'w') as f:
    for atom, vertices in atom_vertices.items():
        formatted_vertices = [f"({v[0]},{v[1]},{v[2]})" for v in vertices if len(v) == 3]
        f.write(' '.join(formatted_vertices) + '\n')

# Step 2: Convert vertices to a single-column format and remove duplicates
with open('1_Al_vertice.txt', 'r') as file:
    unique_elements = set(
        element.strip() for line in file
        for element in line.strip()[1:-1].split(' ')
    )
with open('2_one_col_data.txt', 'w') as file:
    file.write('\n'.join(unique_elements))

# Step 3: Remove brackets and commas from coordinates
with open('2_one_col_data.txt', 'r') as file:
    data = file.read()
modified_data = data.replace('(', '').replace(')', '').replace(',', ' ')
with open('3_modified_one_col_data.txt', 'w') as file:
    file.write(modified_data)

# Step 4: Remove duplicate sites
with open('3_modified_one_col_data.txt', 'r') as file:
    unique_data = set(line.strip() for line in file)
with open('4_unique_data.txt', 'w') as output_file:
    output_file.write('\n'.join(unique_data) + '\n')

# Step 5: Filter out edge and surface sites
atom_positions = np.loadtxt("4_unique_data.txt")
xmin, xmax = xmin_data + 0.5, xmax_data - 0.5
ymin, ymax = ymin_data + 0.5, ymax_data - 0.5
zmin, zmax = zmin_data + 0.5, zmax_data - 0.5
filtered_atoms = atom_positions[
    (atom_positions[:, 0] >= xmin) & (atom_positions[:, 0] <= xmax) &
    (atom_positions[:, 1] >= ymin) & (atom_positions[:, 1] <= ymax) &
    (atom_positions[:, 2] >= zmin) & (atom_positions[:, 2] <= zmax)
]
np.savetxt("5_filtered_atoms.txt", filtered_atoms, fmt="%.6f", delimiter=" ")

# Step 6: Merge close sites
vertices = np.loadtxt('5_filtered_atoms.txt')
clustering = DBSCAN(eps=1.0, min_samples=5).fit(vertices)
merged_vertices = [np.mean(vertices[clustering.labels_ == label], axis=0)
                   for label in set(clustering.labels_) if label != -1]
merged_vertices = np.array(merged_vertices)

# Calculate center and remove closest point to the center
x_min, x_max = np.min(merged_vertices[:, 0]), np.max(merged_vertices[:, 0])
y_min, y_max = np.min(merged_vertices[:, 1]), np.max(merged_vertices[:, 1])
z_min, z_max = np.min(merged_vertices[:, 2]), np.max(merged_vertices[:, 2])
center = np.array([(x_max + x_min) / 2, (y_max + y_min) / 2, (z_max + z_min) / 2])
distances = np.linalg.norm(merged_vertices - center, axis=1)
merged_vertices = np.delete(merged_vertices, np.argmin(distances), axis=0)
np.savetxt('6_merged_sites.txt', merged_vertices, fmt='%.8f')

# Step 7: Add atom ID and type
atom_data = np.loadtxt('Al.txt')
starting_id = atom_data.shape[0] + 1
row_numbers = np.arange(starting_id, starting_id + merged_vertices.shape[0]).reshape(-1, 1)
fixed_number = np.full((merged_vertices.shape[0], 1), 3)
modified_data = np.hstack((row_numbers, fixed_number, merged_vertices))
np.savetxt('7_modified_sites.txt', modified_data, fmt='%d %d %.6f %.6f %.6f')

# Step 8: Modify original data and sort by type
input_filename = '5_16.data'
output_filename = '8_modified_Al.data'
with open(input_filename, 'r') as file:
    header_lines = [next(file) for _ in range(11)]
df = pd.read_csv(input_filename, delim_whitespace=True, header=None, skiprows=11).iloc[:, :8]
df = df.sort_values(by=df.columns[1]).reset_index(drop=True)
df.iloc[:, 0] = range(1, len(df) + 1)
df = df.iloc[:, [0, 1, 2, 3, 4]]
with open(output_filename, 'w') as file:
    file.writelines(header_lines)
    df.to_csv(file, sep='\t', header=False, index=False)

# Step 9: Combine interstitial sites with original data
with open("8_modified_Al.data", "r") as f:
    original_data = f.readlines()
with open("7_modified_sites.txt", "r") as f:
    atom_text_data = f.readlines()
total_atoms = int(atom_text_data[-1].split()[0])
combined_data = original_data + atom_text_data
combined_data[2] = f"{total_atoms} atoms\n"
combined_data[3] = "3 atom types\n"
masses_section = """Masses

1 26.9815  # Al
2 26.9815  # Al
3 58.74  # Ni\n\n"""
combined_data.insert(9, masses_section)
with open("9_combined_structure.data", "w") as f:
    f.writelines(combined_data)

print("Combined data file created as '9_combined_structure.data'.")

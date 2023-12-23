from Bio.PDB import PDBParser
import math

# Atomic radii dictionary
atomic_radii = {"O": 1.52, "P": 1.8, "N": 1.55, "H": 1.2, "C": 1.7}

# Function to calculate the distance between two points
def calculate_distance(coord1, coord2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))

# Count the number of atoms in a chain
def count_atoms(chain):
    return sum(1 for _ in chain.get_atoms())

# Calculate clashes between two residues
def residue_clashes(residue_one, residue_two, threshold):
    clashes = 0
    for atom1 in residue_one.get_atoms():
        for atom2 in residue_two.get_atoms():
            distance = calculate_distance(atom1.coord, atom2.coord)
            if distance <= atomic_radii[atom1.element] + atomic_radii[atom2.element] - threshold:
                clashes += 1
    return clashes

# Calculate the clash score of a structure
def calculate_clash_score(file_path, threshold):
    clash_count = 0
    structure_parser = PDBParser()
    
    structure = structure_parser.get_structure("structure", file_path)
    model = next(structure.get_models())
    chain = next(model.get_chains())

    for i, residue in enumerate(chain):
        for j in range(i + 2, len(chain)):
            clash_count += residue_clashes(chain[i + 1], chain[j + 1], threshold)

    return 1000 * (clash_count / count_atoms(chain))

# Parameters
file_path = "model_4.pdb"
threshold = 0.4

# Calculate and print the clash score
score = calculate_clash_score(file_path, threshold)
print(f"Clash Score: {round(score, 5)}")

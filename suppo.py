from Bio.PDB import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np
import math


def calculate_distance(coordinates1, coordinates2):
    return sum((c1 - c2)**2 for c1, c2 in zip(coordinates1, coordinates2))


def calculate_RMSD(reference_coords, model_coords):
    rmsd = math.sqrt(sum(calculate_distance(ref, mod) for ref, mod in zip(reference_coords, model_coords)) / len(reference_coords))
    return rmsd


def superimpose_coordinates(reference_coords, model_coords):
    superimposer = SVDSuperimposer()
    superimposer.set(reference_coords, model_coords)
    superimposer.run()
    return superimposer.get_transformed()


def get_atom_coordinates(model_residues, reference_residues):
    coords_model = None
    first_iteration = True
    for residue in reference_residues.get_residues():
        for atom in residue.get_atoms():
            atom_coords = model_residues[residue.id[1]][atom.get_id()].get_coord()
            if first_iteration:
                coords_model = np.array(atom_coords)
                first_iteration = False
            else:
                coords_model = np.vstack([coords_model, atom_coords])
    return coords_model


parser = PDBParser()
model_structures = parser.get_structure("Models", "R1107TS081.pdb")
reference_structure = parser.get_structure("Reference", "R1107_reference.pdb")
reference_residues = reference_structure[0]["0"]
reference_coords = get_atom_coordinates(reference_residues, reference_residues)
# Writing the results to a file
with open("results.txt", "w") as result_file:
    for model in model_structures.get_models():
        model_residues = model["0"]
        model_coords = get_atom_coordinates(model_residues, reference_residues)
        superimposed_result = superimpose_coordinates(reference_coords, model_coords)
        rmsd_value = calculate_RMSD(reference_coords, superimposed_result)
        result_file.write(f"Model number {model.id + 1} RMSD: \n{rmsd_value}\n")

"""
Module used for calculating a protein's accessibility to solvent
using the ShrakeRupley method 

"""

# Import necessary modules from Biopython's PDB package
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
import pandas as pd


def shrake_method(name_file, file, dtf):
    # Create a PDBParser object with QUIET mode enabled
    p = PDBParser(QUIET=1)

    # Parse the PDB file to create a structure object
    structure = p.get_structure(name_file, file)

    # Create a ShrakeRupley object for solvent accessibility calculations
    sr = ShrakeRupley()

    # Compute solvent accessibility at the residue level (level="R")
    sr.compute(structure, level="R")

    # Initialize lists to store data
    residue_names = []
    residue_ids = []
    accessible_surfaces = []

    # Iterate through the structure to calculate and store solvent accessibility for each residue
    for model in structure:
        for chain in model:
            for residue in chain:
                # Check if the residue is an amino acid (protein)
                if PDB.is_aa(residue):
                    residue_name = residue.get_resname()  # Get the residue name
                    surface_residue = round(residue.sasa, 2)  # Get the solvent-accessible surface area
                    residue_id = residue.get_id()[1]  # Get the residue ID
                    
                    # Store data in lists
                    residue_names.append(residue_name)
                    residue_ids.append(residue_id)
                    accessible_surfaces.append(surface_residue)

    # Create a Pandas DataFrame from the collected data
    df = dtf.copy()
    df["shrake_accessible_surface"]= accessible_surfaces

    # Compute the total solvent accessibility for the entire protein structure (level="S")
    sr.compute(structure, level="S")
    total_protein_accessibility = round(structure.sasa, 2)

    # Print the total protein accessibility to solvent
    print(f"The protein's total accessibility to solvent is approximately: {total_protein_accessibility} Angstroms")

    # Return the Pandas DataFrame containing residue data
    return df

if __name__ == "__main__":
    import shrake_mth
    print(help(shrake_mth))

"""
The project's main program 

use command in terminal:

    python main.py pdb_name pdb_file.pdb 100

(Replace pdb_name and pdb_file.pdb with the protein's name and file)

"""

import argparse
import atom_coords as ac
import surface_access as sa
import shrake_mth as sm
import plot_comparison as pc
import pandas as pd


if __name__ == "__main__":

    
    PARSER = argparse.ArgumentParser()
    
    PARSER.add_argument("pdb_name",
                       help="The name of the protein",
                       type = str)

    PARSER.add_argument("pdb_file", 
                        help="the pdb file", 
                        type=str)
    
    PARSER.add_argument("nb_points",
                        help="the number of points on the spheres",
                        type=int)
    
    ARGS = PARSER.parse_args()
    
    PDB_NAME = ARGS.pdb_name

    PDB_FILE = ARGS.pdb_file

    NBR_PTS = ARGS.nb_points
    
    print("Step 1 : Reading pdb file and extracting coordinates")

    coa = ac.coords_atoms(PDB_FILE)
    
    print("Step 2 : Calculating accessible surface and total surface for each atom")

    asa = sa.find_accessible_points(coa)

    print("Step 3 : Summing accessible surface of atoms and calculating percentages of accessibility")

    asr = sa.accessibility_surf(asa)

    print(asr)
    
    print("Step 4 : Summing accessible surface of residues")
    
    asp = sa.surf_prot(asr)
    
    print(asp)
	
    print("Step 5 : Calculating accessible surface of residus and protein using the ShrakeRupley method")
    
    sr = sm.shrake_method(PDB_NAME, PDB_FILE, asr)
	
	
    print("Step 6 : Plot of the comparison between our method and the Shrake Method from Bio python")
    pc.comp_plot(sr)

    print("Our plot and dataframe are saved in the /results folder as 'plot.png' and 'access_df.csv':)")
	


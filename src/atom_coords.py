"""
Program that reads a PDB file and extracts each atom coordinates
"""

import pandas as pd


def coords_atoms(file):
    
    """"
    Arguments :
        file : the PDB file
    Return :
        df_coords : Pandas DataFrame of atoms coordinates
    """
    
    # Open the specified file in read mode
    with open(file, "r") as file:
        
        # Create an empty list to store atom coordinates
        list_coords = list()
        
        # Iterate through each line in the file
        for line in file:
            
            # Check if the line starts with "ATOM"
            if line[0:4] == "ATOM":
                # Extract atom name from the line and remove leading/trailing whitespace
                atom_name = line[12:16].strip()
                
                # Check if the atom name does not start with "H" (ignoring hydrogen atoms)
                if not atom_name.startswith("H"):
                    # Create a dictionary to store atom coordinates
                    dict_coords = {}

                    
                    dict_coords["atom"] = str(line[12:16])
                    dict_coords["residue"] = str(line[17:20])
                    dict_coords["residue_nb"] = int(line[22:26])
                    dict_coords["X"] = float(line[30:38])
                    dict_coords["Y"] = float(line[38:46])
                    dict_coords["Z"] = float(line[46:54])

                    # Append the dictionary to the list of coordinates
                    list_coords.append(dict_coords)
                
        # Create a Pandas DataFrame from the list of coordinates
        df_coords = pd.DataFrame(list_coords)
    
    # Return the DataFrame containing atom coordinates
    return df_coords



if __name__ == "__main__":
    # Import the 'coords_atoms' module 
    import atom_coords
    # Display the help/documentation for the 'coords_atoms' module
    print(help(atom_coords))

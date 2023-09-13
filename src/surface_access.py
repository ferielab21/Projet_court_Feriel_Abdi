"""
Module used for selecting neighboors, creating spheres and calculating distances
and surface accessibility

"""

import pandas as pd
import numpy as np


def find_neighbors(index, df, trshld):
    """
    Function that selects neighbors
    
    Arguments:
        index: Line index
        df: Coordinates data frame
        trshld: Distance limit between two neighbors
    
    Returns:
        df_neighbors: Data frame containing the atom's neighbors
    """
    
    # Initialize an empty list to store the neighbors
    list_neighbors = []
    
    # Get the central atom's information from the DataFrame
    atom_central = df.iloc[index]
    
    # Iterate through each line in the DataFrame
    for line in range(len(df)):
        
        # Skip the central atom's row (itself)
        if df.iloc[line].name == index:
            continue
        
        # Get information about a neighboring atom
        atom_n = df.iloc[line]

        # Calculate the distance between the central atom and the neighbor
        distance = np.sqrt((atom_central["X"] - atom_n["X"])**2 +
                            (atom_central["Y"] - atom_n["Y"])**2 +
                            (atom_central["Z"] - atom_n["Z"])**2)

        # Check if the distance is less than or equal to the threshold
        if distance <= trshld:
            # If yes, add the neighbor to the list
            list_neighbors.append(atom_n)

    # Create a DataFrame from the list of neighbors
    df_neighbors = pd.DataFrame(list_neighbors)

    # Return the DataFrame containing the atom's neighbors
    return df_neighbors




def create_sphere(atom, radius, nb_points):
    """
    Function that creates a sphere with n points
    depending on the atom's radius
    
    Arguments:
        atom: The selected atom
        radius: Atom's radius
        nb_points: Number of sphere points
    
    Returns:
        sphere_points: Array with 3D sphere coordinates
    """
    
    # Generate points uniformly distributed on the surface of a sphere
    indexes = np.arange(0, nb_points, dtype=float) + 0.5
    golden_angle = np.pi * (1 + 5**0.5)
    phi = np.arccos(1 - 2 * indexes / nb_points)
    theta = golden_angle * indexes
    
    # Initialize an array to store the 3D coordinates of the points on the sphere
    sphere_points = np.zeros((nb_points, 3))
    
    # Calculate the x, y, and z coordinates of the points on the sphere
    sphere_points[:, 0] = radius * np.cos(theta) * np.sin(phi) + atom["X"]
    sphere_points[:, 1] = radius * np.sin(theta) * np.sin(phi) + atom["Y"]
    sphere_points[:, 2] = radius * np.cos(phi) + atom["Z"]

    # Create a DataFrame to store the coordinates of the points
    df_sphere = pd.DataFrame(sphere_points)

    # Return the array with 3D sphere coordinates
    return sphere_points



def find_accessible_points(dataframe):
    """
    Function that calculates the surface accessibility of atoms
    
    Arguments:
        dataframe: Coordinates data frame
    
    Returns:
        df_surf_access: Data frame containing each atom's surface 
                       accessibility and total surface
    """
    
    # Define van der Waals radii for different atom types
    ATOM_VDW = {
        "C": 2.0,  
        "O": 1.4,  
        "N": 1.5,
        "S": 1.85
    }
    
    # Calculate a threshold distance based on carbon atom's radius and water molecule's radius
    treshold = ATOM_VDW['C'] * 2 + 1.4
    
    # Number of points used to create a sphere
    numb_points = 100
    
    # Lists to store ratios and accessible surface areas
    ratio_tot = []
    accessible_surface_atoms = []  
    surface_atom = []
    
    # Iterate through each atom in the data frame
    for line_index, row in dataframe.iterrows():
        atom_c = row
        atom_name = atom_c["atom"].strip()
        atom_key = atom_name[0]
        
        # Get the van der Waals radius for the atom type and add a constant
        atom_radius = ATOM_VDW.get(atom_key, 0) + 1.4
        
        # Find neighboring atoms within the specified threshold
        neighbors = find_neighbors(line_index, dataframe, treshold)
        
        # Create points on a sphere around the current atom
        sp_points = create_sphere(atom_c, atom_radius, numb_points)

        # List to store solvent-accessible points
        accessible_points = 0

        # Calculate distances between each point on the sphere and each neighbor of the atom
        for coor in sp_points:
            point_X = coor[0]
            point_Y = coor[1]
            point_Z = coor[2]
            is_accessible = True
            for _, row in neighbors.iterrows():
                
                atom_neigh = row
                atom_name_neigh = atom_neigh["atom"].strip()
                atom_key_neigh = atom_name_neigh[0]
                
                # Get the van der Waals radius for the neighbor atom type and add a constant
                radius_neigh = ATOM_VDW.get(atom_key_neigh, 0) + 1.4

                # Calculate the distance between the point and the neighbor atom
                dist = np.sqrt((point_X - atom_neigh["X"])**2 +
                                    (point_Y - atom_neigh["Y"])**2 +
                                    (point_Z - atom_neigh["Z"])**2)

                if dist <= radius_neigh:
                    is_accessible = False
                    break

            if is_accessible:
                accessible_points += 1

        # Calculate the ratio of accessible points to total points on the sphere
        ratio = (accessible_points) / (numb_points)
        ratio_tot.append(ratio)
        
        # Calculate the accessible surface area and the total surface area
        surface = ratio * 4 * np.pi * atom_radius**2
        total_surface = 4 * np.pi * atom_radius**2
        
        accessible_surface_atoms.append(surface)  
        surface_atom.append(total_surface)
    
    # Create a new data frame with added columns for accessible surface and total surface
    df_surf_access = dataframe.copy()
    df_surf_access["Accessible_surface"] = accessible_surface_atoms
    df_surf_access["Total_surface"] = surface_atom
    
    return df_surf_access



def accessibility_surf(df_atom_surface):
    """
    Function that calculates the residues and protein's surface accessibility to solvent
    
    Arguments:
        df_atom_surface: Data frame containing each atom's surface 
                         accessibility and total surface
    
    Returns: 
        df_residue_surface: Data frame containing each residue's surface accessibility
                            and percentage of accessibility
    """
    
    # Create a copy of the input data frame
    residue_access = df_atom_surface.copy()
    
    # Group by 'residue_nb' and 'residue' and sum the other columns
    residue_access = residue_access.groupby(['residue_nb', 'residue']).sum()
    
    # Calculate the 'Percentage_accessibility' column as a percentage of accessible surface area
    residue_access["Percentage_accessibility"] = (residue_access['Accessible_surface'] / residue_access['Total_surface']) * 100
    
    # Create a new data frame with only the desired columns
    df_residue_surface = pd.DataFrame(residue_access, columns=['Accessible_surface', 'Percentage_accessibility'])
    df_residue_surface.to_csv("../results/access_df.csv")
	
    # Return the data frame with residue-level accessibility and the formatted result
    return df_residue_surface


def surf_prot(df_residue):
    """
    Function that sums the total accessible surface of the protein
    
    Arguments : 
        df_residue : Data frame containing each residue's surface accessibility
        result: Protein's surface accessibility
    """
    # Calculate the total accessible surface area of the protein
    df_res_surf = df_residue.copy()
    surface_prot = df_res_surf['Accessible_surface'].sum()
    
    # Format the result as a string
    result = f"The protein's surface accessibility is around {round(surface_prot ,2)} A°^2"
    
    return result


if __name__ == "__main__":
    import surface_access
    print(help(surface_access))
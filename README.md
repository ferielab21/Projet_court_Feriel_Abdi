# Projet_court_Feriel_Abdi
---------------------------------------------------------------
## Description of project :
Calculating the solvent-accessible surface of a protein

## Create conda environment :
Run this command in your terminal to create the conda environment

```bash
$ conda create -n proj_court_env
$ activate proj_court_env
$ conda install pandas
$ conda install -c conda-forge configargparse
$ conda install -c anaconda numpy
$ conda install -c conda-forge biopython
$ conda install -c conda-forge matplotlib-base

# This command will make sure the correct version of biopython is installed
$ pip install biopython --upgrade
```

### Python Packages required :

- Pandas
- BioPython
- Argparse
- Numpy
- Matplotlib

### Data file :

- In the data folder : 1b0q.pdb

### Python scripts : 

In the src folder :
- `main.py` : main program
- `atom_coords.py` : reads a pdb file and extracts the atoms' coordinates
- `surface_access.py` : Selects neighbors, creates spheres, calculates distances and surface accessibility
- `shrake_mth.py` : Calculates the solvent-accessible surface of a protein with ShrakeRupley from BioPython
- `plot_comparison.py` : Compares our method to ShrakeRupley

---
### How to use :

#### Step 1 : Clone the repository (or download the file)
```bash
$ git clone https://github.com/ferielab21/Projet_court_Feriel_Abdi.git
```
### Step 2 : Run the code

```bash
$ python3 main.py pdbname pdbfile nb_of_pts 
```
/!\ Don't forget to replace the arguments with the protein's name, pdb file and the number of points you wish to create the sphere with the correct directories

### Example :
```bash
$ cd Projet_court_Feriel_Abdi
$ python3 src/main.py 1b0q data/1b0q.pdb 96 
```
---
## Author :

Feriel Abdi

M2 - BI
Université Paris Cité

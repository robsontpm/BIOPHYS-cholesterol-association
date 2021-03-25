### where to put output files
# Beware! This should not be the same folder where you store other data
# as usually all data in this folder will be moved to a backup folder
# this can however be a subfolder of your data folder. 
OUTPUT_DIRPATH = "/path/to/output/folder"

### Input file - it needs to be in PDB format unfortunately. Use gromacs to convert your trajectory. ###
PDB_TRAJ_FILEPATH = "/path/to/your/file.pdb"
# TPR is used to run gromacs to compute dipoles
TPR_FILEPATH = "/media/wbbib/bozena/Oksysterole/Gromacs_5/MM_monomer/Chol/md_MM_100ps_prod.tpr"

### definition of the main molecules in the file. We assume they are consecutive atoms. 
# NOTICE for 'plane_atoms_ids':
#   we fit consecutive MOLECULES of this type to this atoms 
#   we want:
#    a) all atoms on the same plane
#    b) line from plane_atoms_ids[1] to plane_atoms_ids[2] on time t-1 be on the line at t )
#    c) position of plane_atoms_ids[2] on t-1 the same as in t)
dipoles_ndx = [2, 3, 5, 8, 9, 11, 14, 16, 18, 19, 23, 26, 29, 32, 35, 36, 40, 42, 45, 48, 50, 52, 56, 59, 62, 65, 67, 71,]
MOLECULES = {'CHL': {'def': {'lo': 1, 'up': 74, 'plane_atoms_ids': [9, 18, 35]}, 'dipoles_ndx': dipoles_ndx}}


# you should define a function of this name
# it should decide if a given atom is a SOLVENT
# or not (it is main molecule then)
def is_solvent(item):
	if isinstance(item, int):
		return (item >= 75)
	else:
		try:
			return item["id"] >= 75
		except:
			return item[0] >= 75


# you should define a function of this name
# it should decide if a given atom is a SOLVENT
# AND a ION. Usefull if we want to distinguish
# those to components.
def is_ion(item):
	if isinstance(item, int):
		return (item >= 6441)
	else:
		try:
			return item["id"] >= 6441
		except:
			return item[0] >= 6441			


# you can use this or write your own function.
def is_water(item):
	return is_solvent(item) and not is_ion(item)


# True == use only heavy atoms for this part of the simulation.
# that is, we discar hydrogen atoms completly from a given molecules
# we have two possible molecules: SOLVENT (so water and ions)
# and MOLECULES, that is, items not marked as SOL in PDB files.
SOLVENT_HAO = True
MOLECULES_HAO = True

# this is just for debug purposes, should be bigger than
# number of your frames in PDB file.
# you can set this to smaller number to exit computation at 
# a given frame.  
MAX_ANALYSIS_FRAMES = 10000000
# this is rarely used for debug purposes. It should not be changed.
MAX_VORO_FRAMES = 0
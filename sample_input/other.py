### where to put output files
# Beware! This should not be the same folder where you store other data
# as usually all data in this folder will be moved to a backup folder
# this can however be a subfolder of your data folder. 
OUTPUT_DIRPATH = "/path/to/output/folder"

### Input file - it needs to be in PDB format unfortunately. Use gromacs to convert your trajectory. ###
PDB_TRAJ_FILEPATH = "/path/to/your/file.pdb"

### definition of the main molecules in the file. We assume they are consecutive atoms. 
# NOTICE for 'plane_atoms_ids':
#   we fit consecutive MOLECULES of this type to this atoms 
#   we want:
#    a) all atoms on the same plane
#    b) line from plane_atoms_ids[1] to plane_atoms_ids[2] on time t-1 be on the line at t )
#    c) position of plane_atoms_ids[2] on t-1 the same as in t)
MOLECULES = {
	'CHL1': {'def': {'lo': 1, 'up': 74, 'plane_atoms_ids': [9, 18, 35]}}, 
	'CHL2': {'def': {'lo': 75, 'up': 148, 'plane_atoms_ids': [83, 92, 109]}}
}
# Used in the 'layers.py' script. 
# how many layers of SOLVENT should be computed streaching out from
# the main molecules. A SOL molecule is at distance 1 if its voronoi cell
# has common face with the molecule itself. Then all molecules that has common
# face with molecules in layer 1 comprise the layer 2, etc.
# for a solvent molecule S and main molecule M we define L[M, S] = number
# of the solvent layer of M in which S is. 
COLLECT_LAYERS = [1,2,3]
# Which cross layers should be computed? 
# Crosslayer 'n' for M1, M2 is defined as all SOL molecules S for which 
# L[M1, S] + L[M2, S] = n
CROSS_LAYERS = [2,3,4]
# Similarly, we compute crosslayers for all molecules.
# basically, you can have any number of main molecules
# but we distinguish only 2 cases: crosslayer of all molecules
# and pairwise crosslayers. 
COMMON_CROSS_LAYERS = [3,4,5,6]
# this is just for debug purposes, should be bigger than
# number of your frames in PDB file.
# you can set this to smaller number to exit computation at 
# a given frame.  
MAX_ANALYSIS_FRAMES = 10000000
# this is rarely used for debug purposes. It should not be changed.
MAX_VORO_FRAMES = 0

# you should define a function of this name
# it should decide if a given atom is a SOLVENT
# or not (it is main molecule then)
def is_solvent(item):
	if isinstance(item, int):
		return (item >= 149)
	else:
		try:
			return item["id"] >= 149
		except:
			return item[0] >= 149


# you should define a function of this name
# it should decide if a given atom is a SOLVENT
# AND a ION. Usefull if we want to distinguish
# those to components.
def is_ion(item):
	if isinstance(item, int):
		return (item >= 6494)
	else:
		try:
			return item["id"] >= 6494
		except:
			return item[0] >= 6494			


# you can use this or write your own function.
def is_water(item):
	return is_solvent(item) and not is_ion(item)


# True == use only heavy atoms for this part of the simulation.
# that is, we discar hydrogen atoms completly from a given molecules
# we have two possible molecules: SOLVENT (so water and ions)
# and MOLECULES, that is, items not marked as SOL in PDB files.
SOLVENT_HAO = True
MOLECULES_HAO = True

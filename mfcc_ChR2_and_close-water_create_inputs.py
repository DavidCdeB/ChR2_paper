import os
import sys
import pyframe

print(f'Using PyFraME version: {pyframe.__version__}')
print(f'Using PyFraME file: {pyframe.__file__}')
print(f'Input file: {sys.argv[1]}')

filename = sys.argv[1]
system = pyframe.MolecularSystem(input_file=filename)

# QM region split:
system.split_fragment_by_identifier(identifier='472_LYR', new_names=['LYS', 'RET'],
fragment_definitions=[['N', 'HN', 'HA', 'C', 'O', 'CA'],
['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3',
 'CE', 'HE2', 'HE3', 'NZ', 'HZ1', 'C1', 'H11', 'C2', 'HC2', 'C3', 'C4',
 'H41', 'H42', 'H43', 'C5', 'H5', 'C6', 'H6', 'C7', 'H7', 'C80', 'C8', 'H81',
 'H82', 'H83', 'C9', 'H9', 'C10', 'H10', 'C11', 'C12', 'C13', 'H131',
 'H132', 'H133', 'C14', 'H141', 'H142', 'C15', 'H151', 'H152',
 'C16', 'H161', 'H162', 'C17', 'C18', 'H181', 'H182', 'H183', 'C19',
 'H191', 'H192', 'H193']])

core = system.get_fragments_by_identifier(identifiers=['472_RET'])
# print (core)
# print (core.values())
# print (core.keys())

# Set the charge=+1 in the Shiff Base:
for atom in core['472_RET'].atoms:
    # print (atom)
    if atom.name == 'NZ':
        atom.charge = 1.0

print ("---------------")
# for atom in core['472_RET'].atoms:
    # print (atom)


# Lipids split:
system.split_fragment_by_name(
    name='POPC',
    new_names=['POCH', 'POCO', 'POCP'],
    fragment_definitions=[['N', 'C13', 'H13A', 'H13B', 'H13C', 'C14', 'H14A', 'H14B', 'H14C', 'C15', 'H15A', 'H15B',
                           'H15C', 'C12', 'H12A', 'H12B', 'C11', 'H11A', 'H11B', 'P', 'O11', 'O12', 'O13', 'O14',
                           'C1', 'HA', 'HB', 'C2', 'HS', 'O21', 'C21', 'O22', 'C3', 'HX', 'HY', 'O31', 'C31', 'O32'
                          ],
                          [
                              'C22', 'H2R', 'H2S', 'C23', 'H3R', 'H3S', 'C24', 'H4R', 'H4S', 'C25', 'H5R', 'H5S',
                              'C26', 'H6R', 'H6S', 'C27', 'H7R', 'H7S', 'C28', 'H8R', 'H8S', 'C29', 'H91', 'C210',
                              'H101', 'C211', 'H11R', 'H11S', 'C212', 'H12R', 'H12S', 'C213', 'H13R', 'H13S',
                              'C214', 'H14R', 'H14S', 'C215', 'H15R', 'H15S', 'C216', 'H16R', 'H16S', 'C217',
                              'H17R', 'H17S', 'C218', 'H18R', 'H18S', 'H18T'
                          ],
                          [
                              'C32', 'H2X', 'H2Y', 'C33', 'H3X', 'H3Y', 'C34', 'H4X', 'H4Y', 'C35', 'H5X', 'H5Y',
                              'C36', 'H6X', 'H6Y', 'C37', 'H7X', 'H7Y', 'C38', 'H8X', 'H8Y', 'C39', 'H9X', 'H9Y',
                              'C310', 'H10X', 'H10Y', 'C311', 'H11X', 'H11Y', 'C312', 'H12X', 'H12Y', 'C313',
                              'H13X', 'H13Y', 'C314', 'H14X', 'H14Y', 'C315', 'H15X', 'H15Y', 'C316', 'H16X',
                              'H16Y', 'H16Z'
                          ]])


lipid = system.get_fragments_by_name(names=['POCH', 'POCO', 'POCP'])
# close_water = system.get_fragments_by_name_and_distance(names=['SOL'], distance=20, reference=core)

# For example, here:
#   ~/Glycine_dimer_Avogadro/Tests/GFP/Branch_mfcc-hb-debug-N-link/10_Feb_Commit_8cf959cd/6_31_plus_G_star/bond_threshold_0.20/hb_water_close_use_mfcc_true/mfcc_hb_potential_hb_water_close.py
# we use:
close_water = system.get_fragments_by_distance_and_name(distance=20.0,
    reference=core,
    use_center_of_mass=False,
    names=['SOL'])

water = system.get_fragments_by_name(names=['SOL'])
ions = system.get_fragments_by_name(names=['SOD', 'CLA'])
ChR2 = system.get_remaining_fragments()

system.set_core_region(fragments=core, basis='6-31+G*')

system.add_region(name='lipid',
                  fragments=lipid,
                  use_standard_potentials=True,
                  standard_potential_model='ALEP',
                  standard_potential_exclusion_type='mfcc',
                  use_mfcc=True,
                  mfcc_order=3)

print(f'Number of lipid: {len(lipid)//3}')

system.add_region(
    name='water', fragments=water,
    use_standard_potentials=True, standard_potential_model='SEP'
)
print(f'Number of water molecules: {len(water)}')

system.add_region(
    name='ions',
    fragments=ions,
    use_standard_potentials=True,
    standard_potential_model="SEP"
)
print(f'Number of ions: {len(ions)}')

system.add_region(
    name='ChR2',
    fragments=ChR2,
    use_mfcc=True,
    use_multipoles=True, # we add this to have folders with 'multipoles' word name.
    multipole_order=2,
    xcfun='CAMB3LYP',
    basis='loprop-6-31+G*',
    use_polarizabilities=True,
    mfcc_order=2
)

# For close-water we don't use use_mfcc=True, since we're not going to use mfcc_hb : 
system.add_region(
    name='close_water',
    fragments=close_water,
    use_multipoles=True,
    multipole_order=2,
    xcfun='CAMB3LYP',
    basis='loprop-6-31+G*',
    use_polarizabilities=True
)
print(f'Number of waters in 20 ang distance from core: {len(close_water)}')
# sys.exit()

project = pyframe.Project()
#export DALTON_TMPDIR=/home/projects/dtu_00024/scratch/$USER
# project.scratch_dir = f'/scratch/{os.environ["SLURM_JOB_ID"]}'
project.scratch_dir = f'/home/projects/dtu_00052/scratch/dcabu'
# project.scratch_dir = f'/home/david/scratch/'
project.jobs_per_node = 1
project.mpi_procs_per_job = 40
project.memory_per_job = 4096 * project.mpi_procs_per_job

project.create_embedding_potential(system, only_create_inputs=True)
# project.create_embedding_potential(system)
# project.write_potential(system)
# project.write_core(system)


from togamgxs.main import TOGA
from togamgxs.groups import Groups


#######################
# Optimization Mode
######################

# If equivalence = True fluxes are included in the ISOXML
# output. This allows MOOSE to calculate basic equivalence factors.
# Please note this feature is NOT compatible with optimization
# at this time since the equivalence fluxes cannot be automatically
# condensed with the rest of the library.

# The variable 'legendre_order' corresponds to the legendre order used in
# the scattering representation.
# If optimize = True, 'legendre_order' corresponds to the max order
# which will be considered in the optimization.

calc_equivalence = False
legendre_order = 1

################# Energy Groups #################

# Some common group structures are specified in the clsas Groups (groups.py)
# Set G.<structure> where <structure> is CASMO, LANL, XMAS, or CUSTOM (no others yet available)
# If CUSTOM please edit the CUSTOM function in groups.py to enter your own groups.

# If optimize = True, all group structures will be used.
# If optimize = False, only the group structure specified in num_groups will be used.
# Select a number of delayed precursor groups, between 1 and 8. Default is 6.

G = Groups()
G.CASMO()
group_structure_info = {"group_structure_library":G,
                        "num_groups":2,
                        "num_delayed_groups":6} # num_delayed_groups can be from 1-8


################# Cross Sections #################

# Of the options listed, include desired reaction types.
# Note: some types will necessitate other types to be included.
# Read the OpenMC errors & warnings during runtime
# to determine which types you may need to add.
# 'absorption'
# 'beta'
# 'chi'
# 'chi-delayed'
# 'chi-prompt'
# 'decay-rate'
# 'delayed-nu-fission'
# 'fission'
# 'inverse-velocity'
# 'kappa-fission'
# 'multiplicity matrix'
# 'nu-fission'
# 'nu-scatter matrix'
# 'prompt-nu-fission'
# 'scatter matrix'
# 'total'
# 'transport'

mgxs_types = ['transport', 'kappa-fission', 'total', 'absorption', 'nu-fission', 'fission',
'scatter matrix', 'nu-scatter matrix', 'multiplicity matrix', 'chi', 'inverse-velocity',
'chi-prompt', 'chi-delayed', 'prompt-nu-fission', 'delayed-nu-fission', 'beta', 'decay-rate']

# Note: Generating cross sections "by_nuclide" is not yet offered by this wrapper.

################ Optimization options #################

# If optimize = True, define group_structures below
# and set iso_histogram and/or angle_histogram, or neither
# and set opt_strategy as 'successive' or 'combinations', defaults to 'combinations'.
# Additionally, opt_tolerance will only choose parameters whose
# final answers are within tolerance of the most accurate run.
# opt_tolerance defaults to 1% of the pcm difference (k_mg - k_ce).

# If optimize = False, the user must set
# quantify_error (True/False): where True makes the resultant mgxs run and
# computes the k-eff bias, plots cross sections and a fission rate error plot.


optimize = False
# If optimize = True, set the following variables:
optimization_options = {'opt_strategy':'successive', # choose from: 'successive', 'combinations'
                        'opt_tolerance':5, # [percent pcm]
                        'iso_histogram':False,
                        'angle_histogram':False}

# If optimize = False, you may use the 'quantify_error' variable:
quantify_error = False


################# Problem Documentation #################

# Any string (including an empty string) is accepted in the following variables:
LibraryName = "Example"
Description = ""
Generator = "Self"
TimeCreated = "2021"

################# MPI options #################

# If running in MPI, choose MPI = True
# Note these mpi args are specially generated for a specific machine.
# If running on a different parallel system, please adapt them to your needs.
# You may need to make modifications to the openmc.run() function in
# mgxs.py and optimize.py

MPI = False
MPI_options = {'machinefile':'mach',
               'numprocs':'16',
               'numthreads':'7',
               'mpicommand':'mpirun'}

# Run TOGA to generate desired MGXS libraries
TOGA(mgxs_types, group_structure_info, legendre_order, optimize, quantify_error,
     calc_equivalence, MPI, LibraryName, Description, Generator, TimeCreated)

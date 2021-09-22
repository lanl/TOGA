# Optimization paramenters: user input required
import re
import sys
import math
import pandas as pd
import openmc
import pdb
import matplotlib.pyplot as plt
from togamgxs.main import TOGA
from togamgxs.groups import Groups
from collections import OrderedDict
import os

sphere_size = sys.argv[1]
num_mesh = sys.argv[2]
num_skip = sys.argv[3]
gen_xs = sys.argv[4].lower() == 'true'

#32grp_alt - 1e8x90 --> 0.99963
starting_groups = [0.0, 9.12e3, 2.48e4, 6.76e4, 8.6517e4, 1.1109e5, 1.42642e5, 1.84e5, 2.35178e5, 3.03e5, 3.87742e5, 4.39369e5, 4.8255e5, 5.64161e5, 6.39279e5, 7.24398e5, 8.23e5, 9.30145e5, 1.05399e6, 1.19433e6, 1.353e6, 1.738e6, 2.232e6, 2.865e6, 3.68e6, 4.72367e6, 6.07e6, 7.79e6, 1.e7, 1.2875e7, 1.65e7, 2.e7]

new_groups = []
estart = 0
#estart = 1e4
eend = 0
#eend = 1.0e7
insert = 2
for ct,e in enumerate(starting_groups):
    if e >= eend:
        new_groups.append(e)
        continue
    if e >= estart:
        new_groups.append(e)
        delta = (starting_groups[ct+1]-e)/(insert+1)
        for i in range(insert): new_groups.append(e+(i+1)*delta)
    else: new_groups.append(e)

if gen_xs:
    print (len(new_groups))
    os.system('sleep 2')

mod = int(num_mesh)%int(num_skip)
if mod != 0: sys.exit('The number of mesh points skipped is not a common denominator of the number of mesh points')
num_shells = int(num_mesh)/int(num_skip)

num_groups = str(len(new_groups)-1)
num_polar = '8'
num_NA = '4'
xs_lib_file = 'mgxs_'+num_mesh+'_'+num_skip+'.xml'

if gen_xs:
    #######################
    # Optimization Mode
    ######################

    # If optimize = True, define GroupStrctures below
    # and set iso_histrogram and/or angle_histogram, or neither
    # and set opt_stratgy as 'successive' or 'combinations', defaults to 'combinations'.
    # Additionally, opt_tolerance will only choose parameters whose
    # final answers are within tolerance of the most accurate run.
    # opt_tolerance defaults to 1% of the pcm difference (k_mg - k_ce).

    # If optimize = False, the user must set
    # 1) quantify_error (True/False): where True makes the resultant mgxs run and
    # computes the k-eff bias, cross section plot, and 2D fission rate error plot.
    # 2) tabulation (OrderedDict): keys correspond to reactor states;
    # reactor states will be generated as all possible
    # permutations of the tabulation keys;
    # NOTE: tabulation is unsed at the moment.

    # The variable 'legendre' corresponds to the legendre order used in
    # the scattering representation.
    # If optimize = True, 'legendre' corresponds to the max order
    # which will be considered in the optimization.

    legendre_order = int(num_NA)
    optimize = False
    # If optimize = True, set the following variables:
    opt_strategy = 'combinations' # choose from: 'successive', 'combinations'
    opt_tolerance = 2 #[percent pcm]
    iso_histogram = False
    angle_histogram = False
    # If optimize = False, you may use the 'quantify_error' variable:
    quantify_error = False
    # Can specify execute=True if OpenMC has not yet been run to generate MGXS,
    # or can set execute=False to use existing OpenMC output files
    execute = True

    # tabulation is required, but it's currently not automated to produce
    # the multiple reactor states like it should in the future.
    tabulation = OrderedDict()
    tabulation['fueltemp'] = [300]

    ################# Energy Groups #################
    # Some common group structures are specified in the clsas Groups (groups.py)
    # Set G.<structure> where <structure> is CASMO or CUSTOM (no others yet available)
    # If CUSTOM please edit the CUSTOM function in groups.py to enter your own groups.

    # If optimize = True, all group structures will be used.
    # If optimize = False, only the group structure specified in num_groups will be used.

    G = Groups()
    #G.LANL()
    G.CUSTOM(new_groups)
    num_groups = len(new_groups)-1

    group_structure_info = {"group_structure_library":G,
                        "num_groups":num_groups,
                        "num_delayed_groups":6} # num_delayed_groups can be from 1-8

    calc_equivalence = False

    #num_groups = 30

    ################# Cross Sections #################
    # Of the options listed, include desired reaction types.
    # Note: some types will necessitate other types to be included.
    # Read the OpenMC errors/warnings during runtime
    # to determine which types you may need to add.
    # 'kappa-fission'
    # 'total'
    # 'transport'
    # 'absorption'
    # 'nu-fission'
    # 'fission'
    # 'scatter matrix'
    # 'nu-scatter matrix'
    # 'multiplicity matrix'
    # 'chi'
    # 'inverse-velocity'
    # 'chi-prompt'
    # 'chi-delayed'
    # 'prompt-nu-fission'
    # 'delayed-nu-fission'
    # 'beta'
    # 'decay-rate'

    mgxs_types = ['transport', 'kappa-fission', 'total', 'absorption', 'nu-fission', 'fission',
    'scatter matrix', 'nu-scatter matrix', 'multiplicity matrix', 'chi', 'inverse-velocity',
    'chi-prompt', 'chi-delayed', 'prompt-nu-fission', 'delayed-nu-fission', 'beta', 'decay-rate']

    # Generating cross sections "by_nuclide" is not yet offered

    ################# Problem Documentation #################
    LibraryName = "Example"
    generator = "ExampleGenerator"
    time_created = "today"

    description = "cross sections for material: "

    ################# MPI options #################
    # If running in MPI, set MPI equal to the dict
    # of MPI options; otherwise, set MPI = False
    # Note these mpi args are specifically for TEBOW
    # If running on a different parallel system,
    # please adapt them to your needs.
    MPI = False
    # MPI = {'machinefile':'mach',
    # 'numprocs':'16',
    # 'numthreads':'7',
    # 'mpirun':'mpirun'}


    ################# Run MGXS #################
    TOGA(mgxs_types, group_structure_info, legendre_order, execute, optimize, quantify_error,
     calc_equivalence, MPI, LibraryName, description, generator, time_created)

    os.system('mv mgxs.xml '+xs_lib_file)

#Make Griffin input file
griffin_out = open('sphere_'+sphere_size+'_'+num_mesh+'.i','w')

block_strng = ' '.join(str(i) for i in range(1,int(num_mesh)+1)).rstrip()
matrl_strng = ''
for i in range(1,int(num_shells)+1):
    for j in range(int(num_skip)): matrl_strng += str(i)+' '

matrl_strng = matrl_strng.rstrip()

griffin_out.write('[Mesh]\n')
griffin_out.write('  [gmg]\n')
griffin_out.write('    type = GeneratedIDMeshGenerator\n')
griffin_out.write('    dim = 1\n')
griffin_out.write('    xmin=0\n')
griffin_out.write('    xmax = '+sphere_size+'\n')
griffin_out.write('    nx = '+num_mesh+'\n')
griffin_out.write('    subdomain_id = \''+block_strng+'\'\n')
griffin_out.write('    material_id = \''+matrl_strng+'\'\n')
griffin_out.write('  []\n')
griffin_out.write('[]\n')
griffin_out.write('\n')
griffin_out.write('[Problem]\n')
griffin_out.write('  coord_type = RSPHERICAL\n')
griffin_out.write('[]\n')
griffin_out.write('\n')
griffin_out.write('[AuxVariables]\n')
griffin_out.write('  [./r_disp]\n')
griffin_out.write('  [../]\n')
griffin_out.write('  [./density_aux]\n')
griffin_out.write('    family = MONOMIAL\n')
griffin_out.write('    order = CONSTANT\n')
griffin_out.write('  [../]\n')
griffin_out.write('  [./temp]\n')
griffin_out.write('    initial_condition = 300\n')
griffin_out.write('  [../]\n')
griffin_out.write('[]\n')
griffin_out.write('\n')
griffin_out.write('[AuxKernels]\n')
griffin_out.write('  [./r_disp]\n')
griffin_out.write('    type = FunctionAux\n')
griffin_out.write('    variable = r_disp\n')
griffin_out.write('    function = r_disp_func\n')
griffin_out.write('  [../]\n')

for i in range(1,int(num_mesh)+1):
    griffin_out.write('  [./density_aux'+str(i)+']\n')
    griffin_out.write('    type = MaterialRealAux\n')
    griffin_out.write('    variable = density_aux\n')
    griffin_out.write('    property = fuel'+str(i)+'_density\n')
    griffin_out.write('    block = '+str(i)+'\n')
    griffin_out.write('  [../]\n')
griffin_out.write('[]\n')

griffin_out.write('[Functions]\n')
griffin_out.write('  [./r_disp_func]\n')
griffin_out.write('    type = ParsedFunction\n')
griffin_out.write('    value = 0\n')
griffin_out.write('  [../]\n')
griffin_out.write('[]\n')
griffin_out.write('\n')
griffin_out.write('[TransportSystems]\n')
griffin_out.write('  particle = neutron\n')
griffin_out.write('  equation_type = eigenvalue\n')
griffin_out.write('\n')
griffin_out.write('  G = '+str(num_groups)+'\n')
griffin_out.write('\n')
griffin_out.write('  VacuumBoundary = \'right\'\n')
griffin_out.write('\n')
griffin_out.write('  [./sn]\n')
griffin_out.write('    scheme = SAAF-CFEM-SN\n')
griffin_out.write('    AQtype = Gauss-Chebyshev\n')
griffin_out.write('    NPolar = '+num_polar+'\n')
griffin_out.write('    NA = '+num_NA+'\n')
griffin_out.write('    fixed_jacobian = true\n')
griffin_out.write('    fission_source_as_material = true\n')
griffin_out.write('    use_displaced_mesh = true\n')
griffin_out.write('    angle_derivative_scheme = SNSweep\n')
griffin_out.write('    verbose = 1\n')
griffin_out.write('    n_delay_groups = 6\n')
griffin_out.write('  [../]\n')
griffin_out.write('[]\n')
griffin_out.write('\n')
griffin_out.write('[Materials]\n')

xs_ct = 0
for i in range(1,int(num_mesh)+1):
    if (i-1)%int(num_skip) == 0: xs_ct += 1
    griffin_out.write('  [./fuel'+str(i)+']\n')
    griffin_out.write('    type = CoupledFeedbackMatIDNeutronicsMaterial\n')
    griffin_out.write('    block = '+str(i)+'\n')
    griffin_out.write('    isotopes = \'pseudo\'\n')
    #griffin_out.write('    fromFile = true\n')
    griffin_out.write('    densities = 1.0\n')
    griffin_out.write('    library_file = '+xs_lib_file+'\n')
    griffin_out.write('    library_name = \'Example\'\n')
    griffin_out.write('    grid_names = \'states\'\n')
    griffin_out.write('    grid_variables = temp\n')
    #griffin_out.write('    material_id = '+str(xs_ct)+'\n')
    griffin_out.write('    displacements = r_disp\n')
    griffin_out.write('  [../]\n')
griffin_out.write('[]\n')
griffin_out.write('\n')
griffin_out.write('[Executioner]\n')
griffin_out.write('  type = NonlinearEigen\n')
#griffin_out.write('  petsc_options_iname = \'-pc_type -pc_hypre_type -ksp_gmres_restart \'\n')
#griffin_out.write('  petsc_options_value = \'hypre boomeramg 100\'\n')
griffin_out.write('  free_power_iterations = 2\n')
griffin_out.write('  nl_abs_tol = 1e-10\n')
griffin_out.write('  l_tol = 1e-4\n')
griffin_out.write('  l_max_its = 50\n')
griffin_out.write('  output_before_normalization = false\n')
griffin_out.write('  output_after_power_iterations = false\n')
griffin_out.write('[]\n')
griffin_out.write('\n')
griffin_out.write('[Outputs]\n')
griffin_out.write('  exodus = true\n')
griffin_out.write('  csv = true\n')
griffin_out.write('  perf_graph = true\n')
griffin_out.write('[]\n')
griffin_out.close()

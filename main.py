import matplotlib.pyplot as plt
import pyrankvote
from pyrankvote import Candidate, Ballot

import openmc
from groups import *
from getfilename import *
from mgxs import *

########## Optimization Mode #######################

# If equivalence = True fluxes are included in the ISOXML
# output. This allows MOOSE to calculate basic equivalence factors.
# Please note this feature is NOT compatible with optimization
# at this time since the equivalence fluxes cannot be automatically
# condensed with the rest of the library.

# If optimize = True, define GroupStructures below
# and set iso_histrogram and/or angle_histogram, or neither
# and set opt_stratgy as 'successive' or 'combinations', defaults to 'combinations'.
# Additionally, opt_tolerance will only choose parameters whose
# final answers are within tolerance of the most accurate run.
# opt_tolerance defaults to 1% of the pcm difference (k_mg - k_ce).

# If optimize = False, the user must set
# quantify_error (True/False): where True makes the resultant mgxs run and
# computes the k-eff bias, plots cross sections and a fission rate error plot.

# The variable 'legendre' corresponds to the legendre order used in
# the scattering representation.
# If optimize = True, 'legendre' corresponds to the max order
# which will be considered in the optimization.

equivalence = False
legendre = 5
optimize = False
# If optimize = True, set the following variables:
opt_strategy = 'successive' # choose from: 'successive', 'combinations'
opt_tolerance = 5 #[percent pcm]
iso_histogram = False
angle_histogram = False
# If optimize = False, you may use the 'quantify_error' variable:
quantify_error = False

################# Energy Groups #################

# Some common group structures are specified in the clsas Groups (groups.py)
# Set G.<structure> where <structure> is CASMO, LANL, XMAS, or CUSTOM (no others yet available)
# If CUSTOM please edit the CUSTOM function in groups.py to enter your own groups.

# If optimize = True, all group structures will be used.
# If optimize = False, only the group structure specified in num_groups will be used.
# Select a number of delayed precursor groups, between 1 and 8. Default is 6.

G = Groups()
G.CASMO()
num_groups = 70
num_delayed_groups = 6 # options: 1-8

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

mgxs_types = ['transport', 'total', 'absorption', 'nu-fission', 'fission',
'scatter matrix', 'nu-scatter matrix', 'multiplicity matrix', 'chi', 'inverse-velocity',
'chi-prompt', 'chi-delayed', 'prompt-nu-fission', 'delayed-nu-fission', 'beta', 'decay-rate']

# Note: Generating cross sections "by_nuclide" is not yet offered by this wrapper.

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
# mgxs.py and optimize.py.

MPI = False
machinefile = 'mach'
numprocs = '16'
numthreads = '7'
mpicommand = 'mpirun'

################# Run MGXS ################# 

XS = MGXS()
fn = FileName()

# tabulation is defined as the list of reactor states ran.
# Please include as many input files as there are reactor states.
# Under these settings, any file ending in "inp.py" is considered an input file,
# however, the user may change this setting in the following line:

fn.get_filename(ending="inp.py")
tabulation = OrderedDict()
tabulation['states'] = list(range(0, len(fn.listfn)))
mgxs_libs = []
n =  0

for f in fn.listfn: 
	f = f.replace('.py', '')
	print(f)
	module = __import__(f)

	materials = module.materials
	u_root = module.u_root
	geom = module.geom
	settings = module.settings
	mesh = module.mesh
	domain_type = module.domain_type
	domain = module.domain
	subdomain_to_plot = module.subdomain_to_plot

	if domain_type == 'material':
		description = "cross sections for material: "
	elif domain_type == 'cell':
		description = "cross sections for cell: "
	elif domain_type == 'mesh':
		description = "cross sections for mesh cell number: "

    #Launching the mgxs module to compute & optimize multi-group cross sections
	if optimize == True and MPI == True:
		mgxs_libs.append(XS.compute_mgxs(n, LibraryName, u_root,
			settings, geom, mesh, domain_type, domain, subdomain_to_plot, mgxs_types,
			G.GroupStructures, num_delayed_groups=num_delayed_groups, optimize=optimize,
			opt_strategy=opt_strategy, opt_tolerance=opt_tolerance, legendre=legendre,
			iso_hist=iso_histogram, angle_hist=angle_histogram, tabulation=tabulation,
			MPI=MPI, machinefile=machinefile, numprocs=numprocs, numthreads=numthreads,
			mpicommand=mpicommand, Description=Description, Generator=Generator,
			TimeCreated=TimeCreated))
	elif optimize == True and MPI == False:
		mgxs_libs.append(XS.compute_mgxs(n, LibraryName, u_root,
			settings, geom, mesh, domain_type, domain, subdomain_to_plot, mgxs_types,
			G.GroupStructures, num_delayed_groups=num_delayed_groups, optimize=optimize,
			opt_strategy=opt_strategy, opt_tolerance=opt_tolerance, legendre=legendre,
			iso_hist=iso_histogram, angle_hist=angle_histogram, tabulation=tabulation,
			Description=Description, Generator=Generator, TimeCreated=TimeCreated))
	elif optimize == False and MPI == True:
		mgxs_libs.append(XS.compute_mgxs(n, LibraryName, u_root,
			settings, geom, mesh, domain_type, domain, subdomain_to_plot, mgxs_types,
			G.GroupStructures, num_groups=num_groups, num_delayed_groups=num_delayed_groups,
			legendre=legendre, quantify_error=quantify_error, tabulation=tabulation,
			MPI=MPI, machinefile=machinefile, numprocs=numprocs, numthreads=numthreads,
			mpicommand=mpicommand, Description=Description, Generator=Generator,
			TimeCreated=TimeCreated))
	elif optimize == False and MPI == False:
		mgxs_libs.append(XS.compute_mgxs(n, LibraryName, u_root,
			settings, geom, mesh, domain_type, domain, subdomain_to_plot, mgxs_types,
			G.GroupStructures, num_groups=num_groups, num_delayed_groups=num_delayed_groups,
			legendre=legendre, quantify_error=quantify_error, tabulation=tabulation,
			Description=Description, Generator=Generator, TimeCreated=TimeCreated))
	n += 1

################# Select XS options #################
# Now that the optimization has occured for each input file separately,
# we optimize over all the best options globally.
print("Selecting best options over all the reactor states...")

if optimize == True and len(XS.choices) > 1:
	options = OrderedDict()
	candidates = []
	ballots = []

	i = 0
	for state in XS.choices:
		for vote in state:
			if vote not in options.values():
				options[str(i)] = vote
			i += 1

	for key in options.keys():
		candidates.append(Candidate(key))
	for state in XS.choices:
		state_ballot = []
		for vote in state:
			for key, value in options.items():
				if vote == value and Candidate(key) not in state_ballot:
					state_ballot.append(Candidate(key))
		ballots.append(Ballot(ranked_candidates=state_ballot))

	election_results = pyrankvote.instant_runoff_voting(candidates,ballots)
	winners = election_results.get_winners()
	selection = options[str(winners[0])]
	print("The best Group/Legendre combination is: ", selection)

elif optimize == True and len(XS.choices) == 1:
	selection = (XS.choices[0][0][0], XS.choices[0][0][1])
	print("The best Group/Legendre combination is: ", selection)
elif optimize == False:
	selection = (num_groups, legendre)
	print("The best Group/Legendre combination is: ", selection)

################# Print to ISOXML #################
# Now that the library has be optimally condensed,
# the cross sections will be re-formatted to an ISOXML
# file readable by MOOSE
print("Condensing libraries and reformatting to ISOXML...")

output = YAKXS(LibraryName, str(selection[0]), equivalence)
new_structure = openmc.mgxs.EnergyGroups(G.GroupStructures[selection[0]])

for i in range(n):
	xsdata_names = []
	for d in mgxs_libs[i].domains:
		xsdata_names.append(d.name)
	mgxs_libs[i] = mgxs_libs[i].get_condensed_library(new_structure)
	if domain_type == 'mesh':
		mgxs_file, new_materials, new_geometry = mgxs_libs[i].create_mg_mode()
	else:
		mgxs_file, new_materials, new_geometry = mgxs_libs[i].create_mg_mode(
			xsdata_names=xsdata_names)
	mgxs_file.name = str(i)+'_mgxs.h5'
	mgxs_file.export_to_hdf5(filename=mgxs_file.name)

output.export_to_xml(selection[0], num_delayed_groups, mgxs_types, selection[1], equivalence,
	mgxs_libs, tabulation, Description=Description, Generator=Generator,
	TimeCreated=TimeCreated)

print("Done reformatting... See mgxs.xml for output.")

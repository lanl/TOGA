# Optimization paramenters: user input required

from collections import OrderedDict

import openmc
import pyrankvote
from pyrankvote import Candidate, Ballot
import matplotlib.pyplot as plt
import pdb

from togamgxs.getfilename import FileName
from togamgxs.mgxs import MGXS
from togamgxs.yakxs import YAKXS



def TOGA(mgxs_types, group_structure_info, legendre_order, execute=True, optimize=False,
         quantify_error=False, calc_equivalence=False, MPI=False,
         LibraryName="MGXS_library", Description="", Generator="Self", TimeCreated=None):

    equivalence = calc_equivalence
    legendre = legendre_order

    if optimize:
        opt_strategy = optimize['opt_strategy']
        opt_tolerance = optimize['opt_tolerance']
        iso_histogram = optimize['iso_histogram']
        angle_histogram = optimize['angle_histogram']
        quantify_error = False

    if MPI:
        machinefile = MPI['machinefile']
        numprocs = MPI['numprocs']
        numthreads = MPI['numthreads']
        mpicommand = MPI['mpirun']

    G = group_structure_info['group_structure_library']
    num_groups = group_structure_info['num_groups']
    num_delayed_groups = group_structure_info['num_delayed_groups']

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
        try:
            mesh = module.mesh
        except AttributeError:
            mesh = None
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
        if optimize == True and MPI:
            mgxs_libs.append(XS.compute_mgxs(n, LibraryName, u_root, settings, geom, mesh, domain_type, domain, subdomain_to_plot, mgxs_types,
                G.group_structures, num_delayed_groups=num_delayed_groups, execute=execute, optimize=optimize, opt_strategy=opt_strategy, opt_tolerance=opt_tolerance, legendre=legendre,
                iso_hist=iso_histogram, angle_hist=angle_histogram, tabulation=tabulation,
                MPI=MPI, machinefile=machinefile, numprocs=numprocs, numthreads=numthreads, mpicommand=mpicommand,
                Description=Description, Generator=Generator, TimeCreated=TimeCreated))
        elif optimize == True and not MPI:
            mgxs_libs.append(XS.compute_mgxs(n, LibraryName, u_root, settings, geom, mesh, domain_type, domain, subdomain_to_plot, mgxs_types,
                G.group_structures, num_delayed_groups=num_delayed_groups, execute=execute, optimize=optimize, opt_strategy=opt_strategy, opt_tolerance=opt_tolerance, legendre=legendre,
                iso_hist=iso_histogram, angle_hist=angle_histogram, tabulation=tabulation,
                Description=Description, Generator=Generator, TimeCreated=TimeCreated))
        elif optimize == False and MPI:
            mgxs_libs.append(XS.compute_mgxs(n, LibraryName, u_root, settings, geom, mesh, domain_type, domain, subdomain_to_plot, mgxs_types,
                G.group_structures, num_groups=num_groups, num_delayed_groups=num_delayed_groups, execute=execute, legendre=legendre, quantify_error=quantify_error, tabulation=tabulation,
                MPI=MPI, machinefile=machinefile, numprocs=numprocs, numthreads=numthreads, mpicommand=mpicommand,
                Description=Description, Generator=Generator, TimeCreated=TimeCreated))
        elif optimize == False and not MPI:
            mgxs_libs.append(XS.compute_mgxs(n, LibraryName, u_root, settings, geom, mesh, domain_type, domain, subdomain_to_plot, mgxs_types,
                G.group_structures, num_groups=num_groups, num_delayed_groups=num_delayed_groups, execute=execute, legendre=legendre, quantify_error=quantify_error, tabulation=tabulation,
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
    new_structure = openmc.mgxs.EnergyGroups(G.group_structures[selection[0]])

    for i in range(n):
        xsdata_names = []
        for d in mgxs_libs[i].domains:
            xsdata_names.append(d.name)
        mgxs_libs[i] = mgxs_libs[i].get_condensed_library(new_structure)
        if domain_type == 'mesh':
            mgxs_file, new_materials, new_geometry = mgxs_libs[i].create_mg_mode()
        else:
            mgxs_file, new_materials, new_geometry = mgxs_libs[i].create_mg_mode(xsdata_names=xsdata_names)
        mgxs_file.name = str(i)+'_mgxs.h5'
        mgxs_file.export_to_hdf5(filename=mgxs_file.name)

    output.export_to_xml(selection[0], num_delayed_groups, mgxs_types, selection[1], equivalence,
        mgxs_libs, tabulation, Description=Description, Generator=Generator,
        TimeCreated=TimeCreated)

    print("Done reformatting... See mgxs.xml for output.")

    ##################################################
    #end class

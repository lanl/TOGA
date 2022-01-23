import copy
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import openmc
from optimize import *
from yakxs import *

class MGXS:

	def __init__(self):

		self.choices = []
		self.biases = []

	def compute_mgxs(self, n, libname, universe, settings, geometry, mesh,
		domain_type, domain, subdomain_to_plot, mgxs_types, GroupStructures,
		num_groups=None, num_delayed_groups=6, optimize=False,
		opt_strategy='successive', opt_tolerance=1.0,
		legendre=0, iso_hist=False, angle_hist=False, quantify_error=False,
		tabulation=None, MPI=False, machinefile='', numprocs='', numthreads='',
		mpicommand='', Description="", Generator='INL', TimeCreated='',
		mgxs_bynuclide=[], mgxs_tablewise =[], mgxs_librarywise =[]):

		print("computing mgxs...")
		os.environ['OPENMC_MG_CROSS_SECTIONS'] = 'mgxs.h5'
		if GroupStructures == {}:
			print("ERROR: GroupStructures requires at least one entry. \
				Define GroupStructures in main.py")

		if optimize == True:
			num_groups = max(GroupStructures.keys())
			group_edges = GroupStructures[num_groups]
		else:
			group_edges = GroupStructures[num_groups]
		energy_groups = openmc.mgxs.EnergyGroups(group_edges)

		tallies_file = openmc.Tallies()
		if quantify_error == True:
			t1 = openmc.Tally(tally_id=100)
			t1.filters.append(openmc.EnergyFilter(group_edges))
			if domain_type == 'cell':
				t1.filters.append(openmc.CellFilter(subdomain_to_plot))
			elif domain_type == 'material':
				t1.filters.append(openmc.MaterialFilter(subdomain_to_plot))
			t1.scores = ['flux']
			tallies_file.append(t1, False)
			t2 = openmc.Tally(tally_id=200)
			t2.filters.append(openmc.MeshFilter(mesh))
			t2.scores = ['fission']
			tallies_file.append(t2, True)
			t3 = openmc.Tally(tally_id=300)

		for i, d in enumerate(domain):
			t = openmc.Tally(tally_id=301+i, name=d.name)
			t.filters.append(openmc.EnergyFilter(group_edges))
			if domain_type == 'cell':
				t.filters.append(openmc.CellFilter(d))
			elif domain_type == 'material':
				t.filters.append(openmc.MaterialFilter(d))
			elif domain_type == 'mesh':
				t.filters.append(openmc.MeshFilter(d))
			t.scores = ['flux']
			tallies_file.append(t, False)

		# Base library
		mgxs_lib = openmc.mgxs.Library(geometry)
		mgxs_lib.energy_groups = energy_groups
		mgxs_lib.mgxs_types = mgxs_types
		mgxs_lib.num_delayed_groups = num_delayed_groups
		mgxs_lib.domain_type = domain_type
		mgxs_lib.domains = domain
		mgxs_lib.legendre_order = legendre
		mgxs_lib.by_nuclide = False
		mgxs_lib.check_library_for_openmc_mgxs()
		mgxs_lib.build_library()
		mgxs_lib.add_to_tallies_file(tallies_file, merge=True)
		
		# Isotropic histogram library
		if iso_hist == True:
			iso_mgxs_lib = openmc.mgxs.Library(geometry)
			iso_mgxs_lib.energy_groups = energy_groups
			iso_mgxs_lib.mgxs_types = mgxs_types
			iso_mgxs_lib.num_delayed_groups = num_delayed_groups
			iso_mgxs_lib.domain_type = domain_type
			iso_mgxs_lib.domains = domain
			iso_mgxs_lib.correction = None
			iso_mgxs_lib.scatter_format = 'histogram'
			iso_mgxs_lib.histogram_bins = 10
			iso_mgxs_lib.by_nuclide = False
			iso_mgxs_lib.check_library_for_openmc_mgxs()
			iso_mgxs_lib.build_library()
			iso_mgxs_lib.add_to_tallies_file(tallies_file, merge=True)
		else:
			iso_mgxs_lib = None

		# Angle-dependent histogram library
		if angle_hist == True:
			angle_mgxs_lib = openmc.mgxs.Library(geometry)
			angle_mgxs_lib.energy_groups = energy_groups
			angle_mgxs_lib.mgxs_types = mgxs_types
			angle_mgxs_lib.num_delayed_groups = num_delayed_groups
			angle_mgxs_lib.domain_type = domain_type
			angle_mgxs_lib.domains = domain
			angle_mgxs_lib.correction = None
			angle_mgxs_lib.scatter_format = 'histogram'
			angle_mgxs_lib.histogram_bins = 10
			angle_mgxs_lib.num_azimuthal = 10
			angle_mgxs_lib.by_nuclide = False
			angle_mgxs_lib.check_library_for_openmc_mgxs()
			angle_mgxs_lib.build_library()
			angle_mgxs_lib.add_to_tallies_file(tallies_file, merge=True)
		else:
			angle_mgxs_lib = None

		tallies_file.export_to_xml()
		if MPI == True:
			openmc.run(mpi_args=[mpicommand,'--bind-to','socket','-np',numprocs,
				'-machinefile',machinefile],threads=numthreads,cwd='.')
		else:
			openmc.run()

		ce_spfile = './'+str(n)+'statepoint_ce.h5'
		os.rename('statepoint.'+str(settings.batches)+'.h5',ce_spfile)
		ce_sumfile = './'+str(n)+'summary_ce.h5'
		os.rename('summary.h5',ce_sumfile)
		SPCE = openmc.StatePoint(ce_spfile, autolink=False)
		su = openmc.Summary(ce_sumfile)
		SPCE.link_with_summary(su)
		print("loading data...")
		mgxs_lib.load_from_statepoint(SPCE)

		# Optimization step
		OPT = Optimize(num_groups, legendre)
		if optimize == True:
			tolerance = 1.0+opt_tolerance/100.0
			if opt_strategy == 'successive':
				self.choices.append(OPT.successive(n, MPI, mpicommand, numprocs,
					numthreads, machinefile, settings, SPCE, GroupStructures, mgxs_lib,
					legendre, iso_mgxs_lib, angle_mgxs_lib, tolerance))
			elif opt_strategy == 'combinations':
				self.choices.append(OPT.combinations(n, MPI, mpicommand, numprocs,
					numthreads, machinefile, settings, SPCE, GroupStructures, mgxs_lib,
					legendre, iso_mgxs_lib, angle_mgxs_lib, tolerance))
			else:
				print("Please set opt_strategy to either 'successive' or 'combinations'.")

		# Error Quantification step
		if quantify_error == True:
			k_ce = SPCE.k_combined
			mgxs_file, new_materials, new_geometry = mgxs_lib.create_mg_mode()
			mgxs_file.export_to_hdf5()
			mgxs_file.name = 'mgxs.h5'
			new_materials.export_to_xml()
			new_geometry.export_to_xml()
			settings.energy_mode = 'multi-group'
			settings.temperature = {'method': 'nearest', 'tolerance': 3000}
			settings.export_to_xml()
			tallies_file = openmc.Tallies()
			tallies_file.append(t2, True)
			tallies_file.export_to_xml()
			print("running in multi-group mode")
			if MPI == True:
				openmc.run(mpi_args=[mpicommand,'--bind-to','socket','-np',numprocs,
					'-machinefile',machinefile],threads=numthreads,cwd='.')
			else:
				openmc.run()

			mg_spfile = './statepoint_mg.h5'
			os.rename('statepoint.'+str(settings.batches)+'.h5',mg_spfile)
			mg_sumfile = './summary_mg.h5'
			os.rename('summary.h5',mg_sumfile)
			SPMG = openmc.StatePoint(mg_spfile, autolink=False)
			su = openmc.Summary(mg_sumfile)
			SPMG.link_with_summary(su)

			k_mg = SPMG.k_combined
			bias = 1.0E5 * (k_mg - k_ce)
			print('Continuous-Energy keff = {0:1.6f}'.format(k_ce))
			print('Multi-Group keff = {0:1.6f}'.format(k_mg))
			print('bias [pcm]: {0:1.1f}'.format(bias.nominal_value))

		    # Plotting cross sections
			print("plotting mgxs...")

			flux_ce = SPCE.get_tally(id=100)
			mg_groups_plt = []
			ce_flux_plt = []
			for i in range(num_groups+1):
				mg_groups_plt.append(group_edges[i])
				if i != num_groups:
					ce_flux_plt.append(flux_ce.mean[i][0][0])
					ce_flux_plt.append(flux_ce.mean[i][0][0])
				if i != 0 and i != num_groups:
					mg_groups_plt.append(group_edges[i])

			if domain_type == 'mesh':
				print("Cannot plot CE/MG cross sections for a material when \
					Mesh domain is used, since a mesh does not know what material is in \
					each cell")
			else: 
				for i, domains in enumerate(mgxs_lib.domains):
					if mgxs_lib.domains[i].name == subdomain_to_plot.name:
						dom = mgxs_lib.domains[i]
						position = i
				if domain_type == 'material':
					if subdomain_to_plot.temperature == None:
						fig = openmc.plot_xs(dom, ['total'])
					else:
						fig = openmc.plot_xs(dom, ['total'],
							temperature=subdomain_to_plot.temperature)
				elif domain_type == 'cell':
					if subdomain_to_plot.temperature == None:
						fig = openmc.plot_xs(dom.fill, ['total'])
					else:
						fig = openmc.plot_xs(dom.fill, ['total'],
							temperature=subdomain_to_plot.temperature)
				openmc.plot_xs(new_materials[position], ['total'], plot_CE=False, 
					mg_cross_sections='mgxs.h5', axis=fig.axes[0])
				fig.axes[0].legend(loc=2, prop={'size': 6}).set_visible(True)
				fig.axes[0].set_title("Cross Sections for "+mgxs_lib.domains[position].name
					+' bias [pcm]: {0:1.1f}'.format(bias.nominal_value))
				ax2 = fig.axes[0].twinx()
				ax2.plot(mg_groups_plt, ce_flux_plt, color='red', label='ce flux')
				ax2.set_ylabel('Flux spectrum [particle-cm per source particle]')
				ax2.set_yscale('log')
				ax2.legend(loc=1, prop={'size': 6})
				plt.savefig('xsplot'+str(n))
				plt.tight_layout()
				plt.close()

			# Plotting 2D reaction rates
			if len(mesh.dimension) != 2:
				print("2D error plot can only be produced if 2D mesh is used")
			else:
				print("plotting 2D error plot...")

				mg_fission_tally = SPMG.get_tally(id=200)
				mg_fission_rates = mg_fission_tally.get_values(scores=['fission'])
				mg_fission_rates.shape = mesh.dimension
				mg_fission_rates /= np.mean(mg_fission_rates[mg_fission_rates > 0.])
				mg_fission_rates[mg_fission_rates == 0.] = np.nan
				ce_fission_tally = SPCE.get_tally(id=200)
				ce_fission_rates = ce_fission_tally.get_values(scores=['fission'])
				ce_fission_rates.shape = mesh.dimension
				ce_fission_rates /= np.mean(ce_fission_rates[ce_fission_rates > 0.])
				ce_fission_rates[ce_fission_rates == 0.] = np.nan
				ratios = []
				for i in range(len(ce_fission_rates)):
					ratios.append(np.divide(mg_fission_rates[i]-ce_fission_rates[i],
						ce_fission_rates[i])[:]*100)

				norm = mpl.colors.Normalize(vmin=0, vmax=20.0)
				plt.imshow(ratios, interpolation='none', cmap='jet', norm=norm, origin='lower')
				plt.title('Relative Fission Rate Error (%)')
				plt.colorbar()
				plt.savefig('FissionRates'+str(n))

			# Printing mgxs, flux, std_dev to text
			if domain_type != 'mesh':
				print("outputing error data to text...")

				total_mgxs = mgxs_lib.get_mgxs(subdomain_to_plot, 'total').get_xs(value='mean')
				std_dev_mgxs = mgxs_lib.get_mgxs(subdomain_to_plot, 'total').get_xs(value='std_dev')
				np.savetxt("mgxs_error"+str(n)+".txt", np.c_[total_mgxs, std_dev_mgxs],
					header='Total MGXS for subdomain        std_dev')

				flux = []
				flux_std_dev = []
				for i in range(num_groups):
					flux.append(flux_ce.mean[i][0][0])
					flux_std_dev.append(flux_ce.std_dev[i][0][0])
					np.savetxt("flux_error"+str(n)+".txt", np.c_[flux, flux_std_dev],
						header='Flux for subdomain         std_dev')

		return mgxs_lib

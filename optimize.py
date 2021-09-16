import copy
import matplotlib as mpl 
import matplotlib.pyplot as plt
import openmc
import os

class Optimize:

####################################################################

	def __init__(self, num_groups, legendre):

		self.G = num_groups
		self.L = legendre

####################################################################

	def successive(self, MPI, mpicommand, numprocs, numthreads, machinefile, 
		settings, SPCE, GroupStructures, mgxs_lib,
		legendre, iso_mgxs_lib, angle_mgxs_lib, tolerance):

		print("This is the successive optimization option. An energy", \
			"optimization takes place using L=0 to select the best group structure G.", \
			"Then, a scattering optimization takes place using G to select the best Legendre expansion L.", \
			"Finally, another energy optimization takes place using L to select the best group structure.")
		ranked_choices = []

		BestGroups, BestBias = self.optimize_energy(MPI, mpicommand, numprocs, numthreads, machinefile,
			settings, SPCE, GroupStructures, mgxs_lib, 
			tolerance, legendre=0)

		ranked_choices.append([BestGroups, 0])

		if legendre > 0:
			BestLegendre = self.optimize_scatter(MPI, mpicommand, numprocs, numthreads, machinefile,
				settings, SPCE, legendre, mgxs_lib, 
				iso_mgxs_lib, angle_mgxs_lib, tolerance, GroupStructures, num_groups=BestGroups)

			ranked_choices.append([BestGroups, BestLegendre])

			BestGroups, BestBias = self.optimize_energy(MPI, mpicommand, numprocs, numthreads, machinefile,
				settings, SPCE, GroupStructures, mgxs_lib, 
				tolerance, legendre=BestLegendre)

			ranked_choices.append([BestGroups, BestLegendre])

		else:
			BestLegendre = 0

		print("Best overall Legendre is:", BestLegendre, "for a group structure of:", BestGroups, "groups")

		ranked_choices = ranked_choices[::-1]
		print(ranked_choices)
		self.G = BestGroups 
		self.L = BestLegendre
		self.B = BestBias
		return ranked_choices

####################################################################

	def combinations(self, MPI, mpicommand, numprocs, numthreads, machinefile, 
		settings, SPCE, GroupStructures, mgxs_lib,
		legendre, iso_mgxs_lib, angle_mgxs_lib, tolerance):

		print("This is the combinations optimization option.", \
			"An energy optimization will take place for every Legendre expansion", \
			"specified in main.py.")

		def sort_by_bias(element):
			return element[0]

		Groups = []
		Bias = []
		Costs = []
		temp_choices = []
		ranked_choices = []

		for L in range(legendre+1):
			TempGroups, TempBias = self.optimize_energy(MPI, mpicommand, numprocs, numthreads, machinefile, 
				settings, SPCE, GroupStructures, mgxs_lib, 
				tolerance, legendre=L)
			Groups.append(TempGroups)
			Bias.append(TempBias)
			Costs.append(self.TempCost*(L+1))
			temp_choices.append([TempBias,TempGroups,L])

		BestLegendre = Costs.index(min(Costs))
		BestGroups = Groups[BestLegendre]
		BestBias = Bias[BestLegendre]

		print("Best overall Legendre is:", BestLegendre, "for a group structure of:", BestGroups, "groups")

		temp_choices.sort(key=sort_by_bias)
		for element in temp_choices:
			ranked_choices.append(temp_choices[1::])
		self.G = BestGroups 
		self.L = BestLegendre
		self.B = BestBias
		return ranked_choices

####################################################################

	def optimize_energy(self, MPI, mpicommand, numprocs, numthreads, machinefile, 
		settings, SPCE, GroupStructures, mgxs_lib, 
		tolerance, legendre=0):

		print("optimizing energy...")

		settings.energy_mode = 'multi-group'
		settings.max_order = legendre
		settings.temperature = {'method': 'nearest', 'tolerance': 3000}
		settings.export_to_xml()
		if os.path.isfile('tallies.xml'):
			os.remove('tallies.xml')
		k_ce = SPCE.k_combined

		groups = []
		biases = []
		stdev = []
		costs = []
		for Ncoarse_grps in GroupStructures:
			print("Running with", Ncoarse_grps, "energy groups and Legendre expansion of", legendre)
			new_structure = openmc.mgxs.EnergyGroups(GroupStructures[Ncoarse_grps])
			new_mgxs_file, new_materials, new_geometry = mgxs_lib.get_condensed_library(
				new_structure).create_mg_mode()
			new_mgxs_file.export_to_hdf5(filename='mgxs.h5')
			new_materials.export_to_xml()
			new_geometry.export_to_xml()
			if MPI == True:
				openmc.run(mpi_args=[mpicommand,'--bind-to','socket','-np',numprocs,'-machinefile',machinefile],threads=numthreads,cwd='.')
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
			cost = Ncoarse_grps * bias
			groups.append(Ncoarse_grps)
			biases.append(abs(bias.nominal_value))
			stdev.append(abs(bias.std_dev))
			costs.append(abs(cost.nominal_value))
			print('Continuous-Energy keff = {0:1.6f}'.format(k_ce))
			print('Multi-Group keff = {0:1.6f}'.format(k_mg))
			print('bias [pcm]: {0:1.1f}'.format(bias.nominal_value))

		# Select optimal group structure 
		costs_sorted = sorted(costs)
		i = 0
		choice = costs.index(costs_sorted[i])
		while biases[choice] > biases[-1]*tolerance:
			i = i + 1
			choice = costs.index(costs_sorted[i])
		BestGroupStructure = GroupStructures[groups[choice]]
		print("The least-cost was acheived with ", groups[choice], " groups.")
		self.TempCost = costs_sorted[i]


		plt.errorbar(groups, biases, stdev, ecolor='k', elinewidth=0.5, capsize=2)
		plt.title("Multi-group bias in k-eff wrt CE mode (L=%i)" %(legendre))
		plt.xlabel("Number of groups")
		plt.ylabel("k-eff difference [pcm]")
		plt.savefig("Group_bias"+str(legendre))
		plt.clf()

		return groups[choice], biases[choice]

####################################################################

	def optimize_scatter(self, MPI, mpicommand, numprocs, numthreads, machinefile, 
		settings, SPCE, legendre, mgxs_lib, 
		iso_mgxs_lib, angle_mgxs_lib, tolerance, GroupStructures, num_groups=0):

		print("optimizing scatter representation...")

		settings.energy_mode = 'multi-group'
		settings.temperature = {'method': 'nearest', 'tolerance': 3000}
		if os.path.isfile('tallies.xml'):
			os.remove('tallies.xml')
		k_ce = SPCE.k_combined

		if num_groups == 0:
			num_groups = max(GroupStructures.keys())

		new_structure = openmc.mgxs.EnergyGroups(GroupStructures[num_groups])
		new_mgxs_file, new_materials, new_geometry = mgxs_lib.get_condensed_library(
			new_structure).create_mg_mode()
		new_mgxs_file.export_to_hdf5(filename='mgxs.h5')
		new_materials.export_to_xml()
		new_geometry.export_to_xml()

		scattering = []
		biases = []
		stdev = []
		costs = []
		for i in range(legendre+1):
			print("Running with scattering order", i, "and", num_groups, "energy groups.")
			settings.max_order = i 
			settings.export_to_xml()
			if MPI == True:
				openmc.run(mpi_args=[mpicommand,'--bind-to','socket','-np',numprocs,'-machinefile',machinefile],threads=numthreads,cwd='.')
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
			cost = (i+1) * bias
			scattering.append(i)
			biases.append(abs(bias.nominal_value))
			stdev.append(abs(bias.std_dev))
			costs.append(abs(cost.nominal_value))
			print('Continuous-Energy keff = {0:1.6f}'.format(k_ce))
			print('Legendre order',i,'keff = {0:1.6f}'.format(k_mg))
			print('bias [pcm]: {0:1.1f}'.format(bias.nominal_value))

		plt.errorbar(scattering, biases, stdev, ecolor='k', elinewidth=0.5, capsize=2, 
			color='b', label="Legendre representation")

		if iso_mgxs_lib is not None:
			iso_mgxs_lib.load_from_statepoint(SPCE)
			iso_mgxs_file, iso_materials_file, iso_geometry_file = iso_mgxs_lib.create_mg_mode()
			iso_materials_file.cross_sections = 'mgxs.h5'
			iso_mgxs_file.export_to_hdf5(filename='mgxs.h5')
			iso_materials_file.export_to_xml()
			iso_geometry_file.export_to_xml()
			if MPI == True:
				openmc.run(mpi_args=[mpicommand,'--bind-to','socket','-np',numprocs,'-machinefile',machinefile],threads=numthreads,cwd='.')
			else:
				openmc.run()

			iso_mg_spfile = './statepoint_mg_iso.h5'
			os.rename('statepoint.' + str(settings.batches) + '.h5', iso_mg_spfile)
			iso_mg_sumfile = './summary_mg_iso.h5'
			os.rename('summary.h5', iso_mg_sumfile)
			iso_mgsp = openmc.StatePoint(iso_mg_spfile, autolink=False)
			iso_mgsum = openmc.Summary(iso_mg_sumfile)
			iso_mgsp.link_with_summary(iso_mgsum)

			iso_mg_keff = iso_mgsp.k_combined
			iso_bias = 1.0e5 * (iso_mg_keff - k_ce)
			x = range(legendre+1)
			iso = [abs(iso_bias.nominal_value)] * (legendre+1)
			print('Continuous-Energy keff = {0:1.6f}'.format(k_ce))
			print('Isotropic histogram keff = {0:1.6f}'.format(iso_mg_keff))
			print('bias for iso [pcm]: {0:1.1f}'.format(iso_bias.nominal_value))
			plt.plot(x, iso, color='r', label="Histogram: isotropic (10 bins)")

		if angle_mgxs_lib is not None:
			angle_mgxs_lib.load_from_statepoint(SPCE)
			angle_mgxs_file, angle_materials_file, angle_geometry_file = angle_mgxs_lib.create_mg_mode()
			angle_materials_file.cross_sections = 'mgxs.h5'
			angle_mgxs_file.export_to_hdf5(filename='mgxs.h5')
			angle_materials_file.export_to_xml()
			angle_geometry_file.export_to_xml()
			if MPI == True:
				openmc.run(mpi_args=[mpicommand,'--bind-to','socket','-np',numprocs,'-machinefile',machinefile],threads=numthreads,cwd='.')
			else:
				openmc.run()

			angle_mg_spfile = './statepoint_mg_angle.h5'
			os.rename('statepoint.' + str(settings.batches) + '.h5', angle_mg_spfile)
			angle_mg_sumfile = './summary_mg_angle.h5'
			os.rename('summary.h5', angle_mg_sumfile)
			angle_mgsp = openmc.StatePoint(angle_mg_spfile, autolink=False)
			angle_mgsum = openmc.Summary(angle_mg_sumfile)
			angle_mgsp.link_with_summary(angle_mgsum)

			angle_mg_keff = angle_mgsp.k_combined
			angle_bias = 1.0e5 * (angle_mg_keff - k_ce)
			x = range(legendre+1)
			angle = [abs(angle_bias.nominal_value)] * (legendre+1)
			print('Continuous-Energy keff = {0:1.6f}'.format(k_ce))
			print('Angle histogram keff = {0:1.6f}'.format(angle_mg_keff))
			print('bias for angle [pcm]: {0:1.1f}'.format(angle_bias.nominal_value))
			plt.plot(x, angle, color='g', label="Histogram: angle-dependent (10 bins, 10 azimuthal")

		# Select optimal legendre order 
		costs_sorted = sorted(costs)
		i = 0
		choice = costs.index(costs_sorted[i])
		while biases[choice] > biases[-1]*tolerance:
			i = i + 1
			choice = costs.index(costs_sorted[i])
		print("The least-cost was acheived with a Legendre expansion order of ", choice, ".")

		plt.legend()
		plt.title("Multi-group bias in k-eff wrt CE mode (G=%i)" %(num_groups))
		plt.xlabel("Legendre scattering order")
		plt.ylabel("k-eff difference [pcm]")
		plt.savefig("Scatter_bias"+str(num_groups))
		plt.clf()

		return choice

####################################################################
# end class 

from datetime import date
from collections import OrderedDict
import re
import pdb

import numpy as np
import pandas as pd
from xml.etree import ElementTree as ET
import mcnptools
from openmc._xml import clean_indentation



class MCNP_MGXS:
    def __init__(self, library_name):
        self.yakxs = ET.Element("YakXs")
        self.root = ET.SubElement(self.yakxs, "Multigroup_Cross_Section_Libraries")
        self.root.set("Name", library_name)

        self.reaction_index = OrderedDict([("Absorption",'g'),
                                           ("DNFraction", 'i'),
                                           ("FissionSpectrum", 'g'),
                                           ("DNSpectrum", 'gi'),
                                           ("DNPlambda",'i'),
                                           ("Fission",'g'),
                                           ("NeutronVelocity",'g'),
                                           ("kappaFission",'g'),
                                           ("nuFission", 'p'),
                                           ("Scattering",'pgl'),
                                           #("nuScattering",'pgl'),
                                           ("Total",'g')])

    def format_mctal_raw_output(self, mctal_filename, mgc_tallynum, fns_tallynum, lcs_tallynum):
        # take in mcnp mctal file, calculate same quantities as OpenMC/Serpent

        mgc_options = ['Flux','InverseVelocity','Total',
                       'Capture','Fission',
                       'nuFission', # Note that this can be either total nu-fission or prompt nu-fission, I think depending on the TOTNU card?
                       'DNnuFission', # Delayed neutron (DN) nu-fission
                       'kappaFission',
                       'Absorption','Scatter']

        tfc = mcnptools.MctalTally.TFC
        mctal_file = mcnptools.Mctal(mctal_filename)

        # extract the mgc tally first
        mgc_tally = mctal_file.GetTally(mgc_tallynum)
        mgc_multiplier_bins = mgc_tally.GetMBins()
        mgc_energy_bins = mgc_tally.GetEBins()
        mgc_cell_bins = mgc_tally.GetFBins()

        mcnp_cell_mgc_results = {}
        for cell_idx, cell_bin in enumerate(mgc_tally.GetFBins()):
            mcnp_cell_mgc_results['cell_{:.0f}'.format(cell_bin)] = {}
            for type_idx, mgc_xs_type in enumerate(mgc_options):
                mcnp_cell_mgc_results['cell_{:.0f}'.format(cell_bin)][mgc_xs_type] = []
                for e_idx, energy_bin in enumerate(mgc_energy_bins):
                    arg = [cell_idx, tfc, tfc, tfc, type_idx, tfc, e_idx, tfc]
                    mcnp_cell_mgc_results['cell_{:.0f}'.format(cell_bin)][mgc_xs_type].append(mgc_tally.GetValue(*arg))

        mgc_dict_of_df = {key:pd.DataFrame(value) for key, value in mcnp_cell_mgc_results.items()}
        mgc_df = pd.concat(mgc_dict_of_df, axis=1)
        mgc_df.index = mgc_energy_bins

        # Calculate the neutron velocity from inverse velocity
        #neutron_vel_subset = mgc_df.xs('NeutronVelocity', level=1, axis=1)
        #mgc_df.loc[:,(slice(None),"NeutronVelocity")] = 1./mgc_df.loc[:,(slice(None),"InverseVelocity")]
        mgc_df = mgc_df.join(1.e8/mgc_df.loc[:, (slice(None), "InverseVelocity")].rename(columns={'InverseVelocity':'NeutronVelocity'})).sort_index(axis=1)
        # Convert kappaFission to J from MeV
        mgc_df.loc[:,(slice(None),'kappaFission')] *= 1.6022e-13
        # Reorder the row MultiIndex to have energies from highest to lowest
        mgc_df.sort_index(axis=0, inplace=True, ascending=False)


        # extract the fns tally to get total and delayed fission neutron spectrum
        fns_tally = mctal_file.GetTally(fns_tallynum)
        fns_energy_bins = fns_tally.GetEBins()
        fns_cell_bins = fns_tally.GetFBins()
        fns_delayed_group_bins = fns_tally.GetTBins()
        fns_named_delayed_group_bins = ['prompt']
        for delayed_bin_idx, delayed_bin in enumerate(fns_delayed_group_bins[1:]):
            fns_named_delayed_group_bins.append('delayed_group{}'.format(1+delayed_bin_idx))

        mcnp_cell_fns_results = {}
        for cell_idx, cell_bin in enumerate(fns_cell_bins):
            mcnp_cell_fns_results['cell_{:.0f}'.format(cell_bin)] = {}
            for delayed_bin_idx, delayed_bin_name in enumerate(fns_named_delayed_group_bins):
                mcnp_cell_fns_results['cell_{:.0f}'.format(cell_bin)][delayed_bin_name] = []
                for e_idx, energy_bin in enumerate(fns_energy_bins):
                    arg = [cell_idx, tfc, tfc, tfc, tfc, tfc, e_idx, delayed_bin_idx]
                    mcnp_cell_fns_results['cell_{:.0f}'.format(cell_bin)][delayed_bin_name].append(fns_tally.GetValue(*arg))


        fns_dict_of_df = {key:pd.DataFrame(value) for key, value in mcnp_cell_fns_results.items()}
        fns_df = pd.concat(fns_dict_of_df, axis=1)
        fns_df.index = fns_energy_bins

        # Calculate the total FNS from the prompt+delayed
        fns_df = fns_df.join(pd.concat([fns_df.sum(level=[0],axis=1)],keys=['total'], axis=1).swaplevel(0,1,axis=1)).sort_index(axis=1)
        # Reorder the row MultiIndex to have energies from highest to lowest
        fns_df.sort_index(axis=0, inplace=True, ascending=False)
        # Calculate the fission spectrums as fractions of each bin, not total
        fns_df = (fns_df.loc[:, (slice(None),slice(None))] / fns_df.sum(axis=0)).fillna(0.0)


        # extract the LCS tally and convert to the scattering matrix
        # analyze the legendre coefficients (LCS)
        lcs_tally = mctal_file.GetTally(lcs_tallynum)
        lcs_incident_energy_bins = lcs_tally.GetUBins()
        lcs_outgoing_energy_bins = lcs_tally.GetEBins() # EBins and UBins might actually need to be swapped: UBins might be outgoing? But I think not, this should be correct.
        lcs_legendre_bins = lcs_tally.GetCBins()
        lcs_cell_bins = lcs_tally.GetFBins()

        mcnp_cell_lcs_results = []
        lcs_column_names = []
        lcs_index_names = [[],[],[]]

        for cell_idx, cell_bin in enumerate(lcs_cell_bins):
            lcs_column_names.append('cell_{:.0f}'.format(cell_bin))
            mcnp_cell_lcs_results.append([])
            for inc_e_idx, inc_e_val in enumerate(lcs_incident_energy_bins):
                for out_e_idx, out_energy_bin in enumerate(lcs_outgoing_energy_bins):
                    for legendre_idx, legendre_order_bin in enumerate(lcs_legendre_bins):
                        if cell_idx == 0:
                            lcs_index_names[0].append(inc_e_val)
                            lcs_index_names[1].append(out_energy_bin)
                            lcs_index_names[2].append(legendre_order_bin)
                        arg = [cell_idx, tfc, inc_e_idx, tfc, tfc, legendre_idx, out_e_idx, tfc]
                        mcnp_cell_lcs_results[-1].append(lcs_tally.GetValue(*arg))

        lcs_nameframe = pd.DataFrame(lcs_index_names, index=['incident_energy','outgoing_energy', 'legendre_order']).T
        lcs_index = pd.MultiIndex.from_frame(lcs_nameframe)
        mcnp_cell_lcs_results = np.array(mcnp_cell_lcs_results)
        lcs_df = pd.DataFrame(data=mcnp_cell_lcs_results, index=lcs_column_names, columns=lcs_index).T
        #lcs_df = lcs_df.drop(1e-11, level='incident_energy').drop(1e-11, level='outgoing_energy')

        # multiply by the scattering cross section of the incident group
        for cell_bin in lcs_df.columns:
            for lcs_incident_energy in lcs_incident_energy_bins:
                for lcs_outgoing_energy in lcs_outgoing_energy_bins:
                    lcs_df.loc[(lcs_incident_energy,lcs_outgoing_energy,slice(1,None)),cell_bin] *=  lcs_df.loc[(lcs_incident_energy,lcs_outgoing_energy,0),cell_bin]
                lcs_df.loc[(lcs_incident_energy,slice(None),slice(None)),cell_bin] *=  mgc_df.loc[lcs_incident_energy,(cell_bin,'Scatter')]

        # determine the legendre order:
        self.legendre_order = len(lcs_legendre_bins) - 1
        self.lcs_legendre_bins = lcs_legendre_bins
        # save the number of groups:
        self.num_groups = len(mgc_energy_bins)
        self.root.set("NGroup", str(self.num_groups))
        # save the number of delayed neutron bins
        self.num_delayed_groups = len(fns_delayed_group_bins[1:])
        # save the lcs group boundaries
        self.lcs_incident_energy_bins = lcs_incident_energy_bins
        self.lcs_outgoing_energy_bins = lcs_outgoing_energy_bins
        # Save the extracted data onto this class
        self.mgc_df = mgc_df
        self.lcs_df = lcs_df
        self.fns_df = fns_df



        if lcs_cell_bins != mgc_cell_bins:
            raise Exception('lcs cell bins {} does not equal mgc cell bins {}, please revise your MCNP tallies and regenerate mctal'.format(lcs_cell_bins, mgc_cell_bins))
        self.domain_bins = [int(item) for item in lcs_cell_bins]
        self.domain_name = 'cell'

        self.tabulation = OrderedDict()
        self.tabulation['states'] = list(range(0, 1))

    def format_outfile_raw_output(self, output_filename):

        precursor_data = {'DNFraction':[],'DNPlambda':[]}

        with open(output_filename, 'r') as f:
            while True:
                line_text = f.readline()
                search_result = re.search(r'precursor    beta-eff', line_text)
                if search_result:
                    # read down two lines
                    f.readline()
                    f.readline()
                    # now, take in the beta fraction and the associated lambda
                    # for each precursor group
                    while True:
                        line_text = f.readline()
                        search_result = re.search(r'\s+[0-9]+\s+(\d+\.\d+)\s+\d+\.\d+\s+\d+\.\d+\s+\d+\.\d+\s+(\d+\.\d+)', line_text)
                        if search_result:
                            precursor_data['DNFraction'].append(float(search_result.group(1)))
                            precursor_data['DNPlambda'].append(float(search_result.group(2)))
                        else:
                            break
                    break
                elif line_text == '':
                    raise Exception("Reached end-of-file for output file {}, meaning delayed neutron data could not be found. Make sure KOPTS card in MCNP input is included, and rerun MCNP to generate this data.".format(output_filename))

            # convert to dataframe
            kopts_df = pd.DataFrame(precursor_data) # If want it to be 1-indexed, use the following: , index=range(1,len(precursor_data['beta'])+1)
            self.kopts_df = kopts_df


    def export_to_csv(self, output_filename_base):
        self.mgc_df.to_csv('{}_mgc.csv'.format(output_filename_base))
        self.lcs_df.to_csv('{}_lcs.csv'.format(output_filename_base))
        self.kopts_df.to_csv('{}_kopts.csv'.format(output_filename_base))


    def export_to_xml(self, output_filename, description='', generator='Self', time_created=None):

        if time_created is None:
            time_created = str(date.today().year)


        # loop over segments
        for id_num in self.domain_bins:
            library = ET.SubElement(self.root, "Multigroup_Cross_Section_Library")
            library.set("ID", str(id_num))
            library.set("Description", "".join([description, self.domain_name, str(id_num)]))
            library.set("Ver", str(1.0))
            library.set("Generator", generator)
            library.set("TimeCreated", time_created)

            # if equivalence:
            #     equivalence_data = ET.SubElement(self.equivalence, "EquivalenceData")
            #     equivalence_data.set("ID", str(id_num))

            ################# Tabulation elements ################
            GridCoordNames = ET.SubElement(library, "Tabulation")
            ReferenceGridIndex = ET.SubElement(library, "ReferenceGridIndex")
            GridCoordNames.text = ""
            ReferenceGridIndex.text = ""
            IndexDict = OrderedDict()

            # if equivalence:
            #     GridCoordNames2 = ET.SubElement(equivalence_data, "Tabulation")
            #     ReferenceGridIndex2 = ET.SubElement(equivalence_data, "ReferenceGridIndex")

            for key in self.tabulation:
                IndexDict[key] = {}
                ReferenceGridIndex.text = \
                ReferenceGridIndex.text+str(len(self.tabulation[key]))+" "
                GridCoordNames.text = GridCoordNames.text+key+" "
                variable = ET.SubElement(library, key)
                # if equivalence:
                #     variable2 = ET.SubElement(equivalence_data, key)
                values = ""
                countvalue = 0
                for value in self.tabulation[key]:
                    countvalue = countvalue+1
                    IndexDict[key].update({value: countvalue})
                    values = values+str(value)+" "
                variable.text = values[:-1]
                # if equivalence:
                #     variable2.text = values[:-1]
            ReferenceGridIndex.text = ReferenceGridIndex.text[:-1]
            GridCoordNames.text = GridCoordNames.text[:-1]
            # if equivalence:
            #     GridCoordNames2.text = GridCoordNames.text
            #     ReferenceGridIndex2.text = ReferenceGridIndex.text

            AllReactions = ET.SubElement(library, "AllReactions")
            AllReactions.text = " ".join(self.reaction_index.keys())

            TablewiseReactions = ET.SubElement(library, "TablewiseReactions")
            TablewiseReactions.text = " ".join(self.reaction_index.keys())

            # If want to loop over multiple MCNP output files at different statepoints, can add a loop here
            # But for now, no loop, so set num_files=1
            num_files = 1

            Table = ET.SubElement(library, "Table")
            Table.set("gridIndex", str(num_files))

            Isotope = ET.SubElement(Table, "Isotope")
            Isotope.set("Name", "pseudo")
            Isotope.set("L", str(self.legendre_order))
            Isotope.set("I", str(self.num_delayed_groups))

            Tablewise = ET.SubElement(Table, "Tablewise")
            Tablewise.set("L", str(self.legendre_order))
            Tablewise.set("I", str(self.num_delayed_groups))

            for reaction_type in self.reaction_index:
                Reaction = ET.SubElement(Tablewise, reaction_type)
                Reaction.set("index", self.reaction_index[reaction_type])
                Reaction.text = ""

                # Now, actually take data from the MCNP tally output dataframe
                #  and output it to xml
                if reaction_type == 'Scattering':
                    Reaction.set("profile", str(1))
                    Profile = ET.SubElement(Reaction, "Profile")
                    Value = ET.SubElement(Reaction, "Value")

                    #g_min = list(datafile[name][groups[0]][rxn]['g_min'])
                    #g_max = list(datafile[name][groups[0]][rxn]['g_max'])
                    #data = list(datafile[name][groups[0]][rxn]['scatter_matrix'])

                    Profile.text = "\n"
                    Value.text = "\n"

                    for legendre_order in self.lcs_legendre_bins:
                        for group_idx in range(1, self.num_groups+1):
                            Profile.text += "1 {}\n".format(self.num_groups)
                    for legendre_order in self.lcs_legendre_bins:
                        #count = 0
                        for group_in in self.lcs_incident_energy_bins[::-1]:
                            for group_out in self.lcs_outgoing_energy_bins[::-1]:
                                Value.text += "{} ".format(self.lcs_df.loc[(group_in, group_out, legendre_order),'cell_{}'.format(id_num)]) #Value.text+str(data[j+(legendre_order+1)*count])+" "
                                #count = count+1
                            Value.text += "\n"
                elif reaction_type == 'FissionSpectrum':
                    Reaction.text += re.sub(r'[\[\]]', '',np.array2string(self.fns_df.loc[:,('cell_{}'.format(id_num), 'prompt')].values))
                elif reaction_type == 'DNSpectrum':
                    Reaction.text += re.sub(r'[\[\]\r\n]', '',np.array2string(self.fns_df.loc[:,('cell_{}'.format(id_num), slice(None,'delayed_group{}'.format(self.num_delayed_groups)))].T.values))
                elif reaction_type == 'DNPlambda' or reaction_type == 'DNFraction':
                    Reaction.text += re.sub(r'[\[\]]', '',np.array2string(self.kopts_df.loc[:,reaction_type].values))
                else:
                    Reaction.text += re.sub(r'[\[\]]', '',np.array2string(self.mgc_df.loc[:,('cell_{}'.format(id_num), reaction_type)].values))


            # If wanted to loop over multiple MCNP output files, this is where that loop would end (before the next few lines)
            Librarywise = ET.SubElement(library, "Librarywise")
            Librarywise.set("L", str(self.legendre_order))
            Librarywise.set("I", str(self.num_delayed_groups))

            clean_indentation(self.yakxs)
            tree = ET.ElementTree(self.yakxs)
            tree.write(output_filename, xml_declaration=True, encoding='utf-8')

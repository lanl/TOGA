

from collections import OrderedDict
from itertools import product
import pdb

import h5py
import numpy as np
import openmc
from openmc._xml import clean_indentation
from xml.etree import ElementTree as ET

from togamgxs.getfilename import FileName

class YAKXS:

####################################################################

    def __init__(self, LibraryName, NGroup, equivalence):

        self.yakxs = ET.Element("YakXs")
        self.root = ET.SubElement(self.yakxs, "Multigroup_Cross_Section_Libraries")
        self.root.set("Name", LibraryName)
        self.root.set("NGroup", NGroup)
        if equivalence:
            self.equivalence = ET.SubElement(self.yakxs, "Equivalence_Data_Library")
            self.equivalence.set("Name", LibraryName)
            self.equivalence.set("NGroup", NGroup)

####################################################################

    def export_to_xml(self, NGroups, num_delayed_groups, mgxs_types, legendre_order, equivalence,
        mgxs_libs, tabulation,
        Description="", Generator="", TimeCreated="",
        mgxs_bynuclide=[], mgxs_tablewise =[], mgxs_librarywise =[]):

        hf = h5py.File('0_mgxs.h5', 'r')
        ID = 0

        fn = FileName()
        fn.get_filename(ending="_mgxs.h5")

        if equivalence:
            fn2 = FileName()
            fn3 = FileName()
            fn2.get_filename(ending="statepoint_ce.h5")
            fn3.get_filename(ending="summary_ce.h5")

        for name in list(hf):
            ID += 1

            groups = list(hf[name])
            reactions = list(hf[name][groups[0]])
            mgxs_tablewise = reactions
            #if 'transport' in mgxs_types:
            #    reactions.append('transport')

            ################# Header elements ################
            library = ET.SubElement(self.root, "Multigroup_Cross_Section_Library")
            library.set("ID", str(ID))
            library.set("Description", Description+name)
            library.set("Ver", str(1.0))
            library.set("Generator", Generator)
            library.set("TimeCreated", TimeCreated)

            if equivalence:
                equivalence_data = ET.SubElement(self.equivalence, "EquivalenceData")
                equivalence_data.set("ID", str(ID))

            ################# Tabulation elements ################
            GridCoordNames = ET.SubElement(library, "Tabulation")
            ReferenceGridIndex = ET.SubElement(library, "ReferenceGridIndex")
            GridCoordNames.text = ""
            ReferenceGridIndex.text = ""
            IndexDict = OrderedDict()

            if equivalence:
                GridCoordNames2 = ET.SubElement(equivalence_data, "Tabulation")
                ReferenceGridIndex2 = ET.SubElement(equivalence_data, "ReferenceGridIndex")

            for key in tabulation:
                IndexDict[key] = {}
                ReferenceGridIndex.text = \
                ReferenceGridIndex.text+str(len(tabulation[key]))+" "
                GridCoordNames.text = GridCoordNames.text+key+" "
                variable = ET.SubElement(library, key)
                if equivalence:
                    variable2 = ET.SubElement(equivalence_data, key)
                values = ""
                countvalue = 0
                for value in tabulation[key]:
                    countvalue = countvalue+1
                    IndexDict[key].update({value: countvalue})
                    values = values+str(value)+" "
                variable.text = values[:-1]
                if equivalence:
                    variable2.text = values[:-1]
            ReferenceGridIndex.text = ReferenceGridIndex.text[:-1]
            GridCoordNames.text = GridCoordNames.text[:-1]
            if equivalence:
                GridCoordNames2.text = GridCoordNames.text
                ReferenceGridIndex2.text = ReferenceGridIndex.text

            AllReactions = ET.SubElement(library, "AllReactions")
            AllReactions.text = self.ReactionTypes(reactions)

            TablewiseReactions = ET.SubElement(library, "TablewiseReactions")
            TablewiseReactions.text = self.ReactionTypes(mgxs_tablewise)

            ################# Tables: permutations of tabulation ################
            for i, f in enumerate(fn.listfn):

                Table = ET.SubElement(library, "Table")
                Table.set("gridIndex", str(i+1))
                if equivalence:
                    Table2 = ET.SubElement(equivalence_data, "Table")
                    Table2.set("gridIndex", str(i+1))
                    Flux = ET.SubElement(Table2, "Flux")
                    Flux.text = ""

                    SPCE = openmc.StatePoint(fn2.listfn[i])
                    if mgxs_libs[i].domain_type == 'mesh':
                        t = SPCE.get_tally(id=301)
                        for E in range(NGroups):
                            result = t.mean[ID-1 + E*len(list(hf))]
                            Flux.text = Flux.text+str(result[0][0])+" "
                    else:
                        t = SPCE.get_tally(name=name)
                        for result in t.mean:
                            Flux.text = Flux.text+str(result[0][0])+" "
                    Flux.text = Flux.text[:-1]

                Isotope = ET.SubElement(Table, "Isotope")
                Isotope.set("Name", "pseudo")
                Isotope.set("L", str(legendre_order))
                Isotope.set("I", str(num_delayed_groups))

                Tablewise = ET.SubElement(Table, "Tablewise")
                Tablewise.set("L", str(legendre_order))
                Tablewise.set("I", str(num_delayed_groups))

                datafile = h5py.File(f,'r')
                for rxn in reactions:

                    if rxn == 'transport':
                        # if rxn == 'total':
                        #     Reaction = ET.SubElement(Tablewise, 'Total')
                        # elif rxn == 'transport':
                        Reaction = ET.SubElement(Tablewise, 'Transport')
                        Reaction.set("index", 'g')
                        Reaction.text = ""
                        if mgxs_libs[i].domain_type == 'mesh':
                            data = mgxs_libs[i].get_mgxs(mgxs_libs[i].domains[0], rxn).get_xs()[ID-1]
                        else:
                            for j, domain in enumerate(mgxs_libs[i].domains):
                                if j == ID-1:
                                    data = mgxs_libs[i].get_mgxs(domain, rxn).get_xs()
                        for d in data:
                            Reaction.text = Reaction.text+str(d)+" "
                        Reaction.text = Reaction.text[:-1]

                    if rxn in self.Conversion:
                        Reaction = ET.SubElement(Tablewise, self.Conversion[rxn])
                        Reaction.set("index", self.DefaultIndex[rxn])
                        Reaction.text = ""

                        # For scattering data
                        if rxn == 'scatter_data':
                            Reaction.set("profile", str(1))
                            Profile = ET.SubElement(Reaction, "Profile")
                            Value = ET.SubElement(Reaction, "Value")

                            g_min = list(datafile[name][groups[0]][rxn]['g_min'])
                            g_max = list(datafile[name][groups[0]][rxn]['g_max'])
                            data = list(datafile[name][groups[0]][rxn]['scatter_matrix'])

                            Profile.text = "\n"
                            Value.text = "\n"

                            for j in range(legendre_order+1):
                                for k in range(len(g_min)):
                                    Profile.text = Profile.text+str(g_min[k])+" "+str(g_max[k])+"\n"
                            for j in range(legendre_order+1):
                                count = 0
                                for g_in in range(NGroups):
                                    for g_out in range(g_min[g_in]-1,g_max[g_in]):
                                        Value.text = Value.text+str(data[j+(legendre_order+1)*count])+" "
                                        count = count+1
                                    Value.text = Value.text+"\n"

                        #For non-scattering data
                        else:
                            #pdb.set_trace()
                            data = list(datafile[name][groups[0]][rxn])
                            for d in data:
                                if rxn == 'inverse-velocity':
                                    Reaction.text = Reaction.text+str(1/d)+" "
                                elif rxn == 'chi-delayed' or rxn  == 'beta':
                                    for g in d:
                                        Reaction.text = Reaction.text+str(g)+" "
                                else:
                                    Reaction.text = Reaction.text+str(d)+" "
                        Reaction.text = Reaction.text[:-1]

            Librarywise = ET.SubElement(library, "Librarywise")
            Librarywise.set("L", str(legendre_order))
            Librarywise.set("I", str(num_delayed_groups))

            clean_indentation(self.yakxs)
            tree = ET.ElementTree(self.yakxs)
            tree.write('./mgxs.xml', xml_declaration=True, encoding='utf-8')

####################################################################

    def ReactionTypes(self, mgxs_types):

        self.Conversion = {
        "absorption": 'Absorption',
        "fission": 'Fission',
        "scatter matrix": 'Scattering',
        "scatter_data": 'Scattering',
        "nu-fission": 'nuFission',
        "kappa-fission": 'kappaFission',
        "chi": 'FissionSpectrum',
        "chi-delayed": 'DNSpectrum',
        "inverse-velocity": 'NeutronVelocity',
        "decay-rate": 'DNPlambda',
        "beta": 'DNFraction',
        "total":'Total'
        }

        self.DefaultIndex = {
        "absorption": 'g',
        "fission": 'g',
        "scatter matrix": 'pgl',
        "scatter_data": 'pgl',
        "nu-fission": 'p',
        "kappa-fission": 'g',
        "chi": 'g',
        "chi-delayed": 'gi',
        "inverse-velocity": 'g',
        "decay-rate": 'i',
        "beta": "pi",
        'total':'g'
        }

        text = ""
        for mgxs in mgxs_types:
            if mgxs in self.Conversion:
                text = text+self.Conversion[mgxs]+" "
            # if mgxs == 'total':
            #     text = text+"Total"+" "
            if mgxs == 'transport':
                text = text+"Transport"+" "
        text = text[:-1]

        return text

####################################################################
#end class

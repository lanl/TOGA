# OpenMC input deck for testbed reactor model

################################################
#                     Imports
################################################

import openmc
import matplotlib.pyplot as plt
import openmc.checkvalue as cv


################################################
#                 Define Materials
################################################

fuel_mat = openmc.Material(material_id=1, name='fuel_mat')
fuel_mat.temperature = 295.
fuel_mat.add_nuclide('U235',0.2)
fuel_mat.add_nuclide('U238',0.8)
fuel_mat.add_nuclide('H1',1.0)
fuel_mat.add_element('Zr',1.0)
fuel_mat.set_density('g/cm3',9.8)
# fuel_mat.add_s_alpha_beta('c_Zr_in_ZrH')
# fuel_mat.add_s_alpha_beta('c_H_in_ZrH')

mod_mat = openmc.Material(material_id=2, name='mod_mat')
mod_mat.add_element('Zr',1.0)
mod_mat.add_nuclide('H1',2.0)
mod_mat.set_density('g/cm3',5.6)
# mod_mat.add_s_alpha_beta('c_Zr_in_ZrH')
# mod_mat.add_s_alpha_beta('c_H_in_ZrH')

vapor_mat = openmc.Material(material_id=3, name='vapor_mat')
vapor_mat.add_nuclide('Na23',1.0)
vapor_mat.set_density('g/cm3',0.1)

pipe_mat = openmc.Material(material_id=4, name='pipe_mat')
pipe_mat.add_element('Mo',1.0)
pipe_mat.set_density('g/cm3',8.0)

monolith_mat = openmc.Material(material_id=5, name='monolith_mat')
monolith_mat.add_nuclide('C12',1.0)
monolith_mat.set_density('g/cm3',1.75)
# monolith_mat.add_s_alpha_beta('c_Graphite')

reflector_mat = openmc.Material(material_id=6, name='reflector_mat')
reflector_mat.add_nuclide('Be9',1.0)
reflector_mat.add_nuclide('O16',1.0)
reflector_mat.set_density('g/cm3',2.9)
# reflector_mat.add_s_alpha_beta('c_Be_in_BeO')
# reflector_mat.add_s_alpha_beta('c_O_in_BeO')

materials = openmc.Materials([fuel_mat,mod_mat,vapor_mat,pipe_mat,monolith_mat,reflector_mat])
materials.export_to_xml()

################################################
#                 Define Geometry
################################################

################# Surfaces #################

# ZCylinders to model pins
pipe_inner = openmc.ZCylinder(surface_id=100, x0=0, y0=0, r=0.6)
pipe_outer = openmc.ZCylinder(surface_id=101, x0=0, y0=0, r=0.7)
mod_r      = openmc.ZCylinder(surface_id=102, x0=0, y0=0, r=0.75)
fuel_r     = openmc.ZCylinder(surface_id=103, x0=0, y0=0, r=0.8)

# top & bottom of the assembly
assembly_z0 = openmc.ZPlane(surface_id=300, z0=-1, boundary_type='reflective')
assembly_z1 = openmc.ZPlane(surface_id=301, z0=1, boundary_type='reflective')

# assembly hexagon
assembly = openmc.model.hexagonal_prism(edge_length=16.166, orientation='y')

# reflector hexagon
reflector = openmc.model.hexagonal_prism(edge_length=27.713, orientation='y',
    boundary_type='vacuum')

################# Cells #################

pipe_cell_inner = openmc.Cell(name= 'pipe_cell_inner', cell_id=100)
pipe_cell_wall = openmc.Cell(name= 'pipe_cell_wall', cell_id=101)
pipe_cell_outer = openmc.Cell(name= 'pipe_cell_outer', cell_id=102)
mod_cell_inner = openmc.Cell(name= 'mod_cell_inner', cell_id=103)
mod_cell_outer = openmc.Cell(name= 'mod_cell_outer', cell_id=104)
fuel_cell_inner = openmc.Cell(name= 'fuel_cell_inner', cell_id=105)
fuel_cell_outer = openmc.Cell(name= 'fuel_cell_outer', cell_id=106)

monolith_cell = openmc.Cell(name= 'monolith_cell', cell_id=201)
assembly_cell = openmc.Cell(name= 'assembly_cell', cell_id=301)

reflect_cell = openmc.Cell(name= 'reflect_cell', cell_id=401)

#################  Regions #################

pipe_cell_inner.region = -pipe_inner & -assembly_z1 & +assembly_z0
pipe_cell_wall.region = -pipe_outer & + pipe_inner & -assembly_z1 & +assembly_z0
pipe_cell_outer.region = +pipe_outer & -assembly_z1 & +assembly_z0
mod_cell_inner.region = -mod_r & -assembly_z1 & +assembly_z0
mod_cell_outer.region = +mod_r & -assembly_z1 & +assembly_z0
fuel_cell_inner.region = -fuel_r & -assembly_z1 & +assembly_z0
fuel_cell_outer.region = +fuel_r & -assembly_z1 & +assembly_z0

assembly_cell.region = assembly & -assembly_z1 & +assembly_z0

reflect_cell.region = ~assembly & reflector & -assembly_z1 & +assembly_z0

# fill cells with material
pipe_cell_inner.fill = vapor_mat
pipe_cell_wall.fill = pipe_mat
pipe_cell_outer.fill = monolith_mat
mod_cell_inner.fill = mod_mat
mod_cell_outer.fill = monolith_mat
fuel_cell_inner.fill = fuel_mat
fuel_cell_outer.fill = monolith_mat
monolith_cell.fill = monolith_mat

reflect_cell.fill = reflector_mat

################# Universes #################

u_pipe = openmc.Universe(universe_id=100)
u_pipe.add_cells(
    [pipe_cell_inner, pipe_cell_wall, pipe_cell_outer])

u_mod = openmc.Universe(universe_id=101)
u_mod.add_cells(
    [mod_cell_inner, mod_cell_outer])

u_fuel = openmc.Universe(universe_id=102)
u_fuel.add_cells(
    [fuel_cell_inner, fuel_cell_outer])

u_monolith = openmc.Universe(universe_id=200,
    cells=(monolith_cell,))

u_reflect = openmc.Universe(universe_id=300,
    cells=(reflect_cell,))

u_root = openmc.Universe(universe_id=400)

################# Lattice #################

lattice = openmc.HexLattice()
lattice.center = (0.0,0.0)
lattice.pitch = (1.88,)
lattice.outer = u_monolith

ring1 = [u_pipe]
ring2 = [u_fuel,u_mod,u_fuel,
         u_mod,u_fuel,u_mod]
ring3 = [u_mod,u_pipe,u_fuel,
         u_pipe,u_mod,u_pipe,
         u_fuel,u_pipe,u_mod,
         u_pipe,u_fuel,u_pipe]
ring4 = [u_pipe,u_fuel,u_mod,
         u_pipe,u_mod,u_fuel,
         u_pipe,u_fuel,u_mod,
         u_pipe,u_mod,u_fuel,
         u_pipe,u_fuel,u_mod,
         u_pipe,u_mod,u_fuel]
ring5 = [u_fuel,u_mod,u_pipe,
         u_fuel,u_mod,u_fuel,
         u_pipe,u_mod,u_fuel,
         u_mod,u_pipe,u_fuel,
         u_mod,u_fuel,u_pipe,
         u_mod,u_fuel,u_mod,
         u_pipe,u_fuel,u_mod,
         u_fuel,u_pipe,u_mod]
ring6 = [u_mod,u_pipe,u_fuel,
         u_mod,u_pipe,u_fuel,
         u_pipe,u_mod,u_fuel,
         u_pipe,u_mod,u_pipe,
         u_fuel,u_mod,u_pipe,
         u_fuel,u_pipe,u_mod,
         u_fuel,u_pipe,u_mod,
         u_pipe,u_fuel,u_mod,
         u_pipe,u_fuel,u_pipe,
         u_mod,u_fuel,u_pipe]
ring7 = [u_pipe,u_fuel,u_mod,
         u_pipe,u_fuel,u_mod,
         u_pipe,u_mod,u_fuel,
         u_pipe,u_mod,u_fuel,
         u_pipe,u_fuel,u_mod,
         u_pipe,u_fuel,u_mod,
         u_pipe,u_mod,u_fuel,
         u_pipe,u_mod,u_fuel,
         u_pipe,u_fuel,u_mod,
         u_pipe,u_fuel,u_mod,
         u_pipe,u_mod,u_fuel,
         u_pipe,u_mod,u_fuel]
ring8 = [u_fuel,u_mod,u_pipe,
         u_fuel,u_mod,u_pipe,
         u_fuel,u_mod,u_fuel,
         u_pipe,u_mod,u_fuel,
         u_pipe,u_mod,u_fuel,
         u_mod,u_pipe,u_fuel,
         u_mod,u_pipe,u_fuel,
         u_mod,u_fuel,u_pipe,
         u_mod,u_fuel,u_pipe,
         u_mod,u_fuel,u_mod,
         u_pipe,u_fuel,u_mod,
         u_pipe,u_fuel,u_mod,
         u_fuel,u_pipe,u_mod,
         u_fuel,u_pipe,u_mod]
ring9=[u_mod]*48

lattice.universes = \
[ring9,ring8,ring7,ring6,ring5,ring4,ring3,ring2,ring1]
assembly_cell.fill = lattice
u_root.add_cells(
    [assembly_cell, reflect_cell])

print(lattice)
geom = openmc.Geometry(u_root)
geom.export_to_xml()

########################
# Plot
########################

# xy plot
# plt.figure(figsize=(12,12))
# u_root.plot(origin=(0,0,0),width=(50,60),color_by='material',pixels=[400,400])
# plt.savefig('Plot')
# plt.clf()

#######################
# Run Settings
#######################

uniform_dist = openmc.stats.Box([-24,-24,-1],[24,24,1],only_fissionable=True)

settings = openmc.Settings()
settings.output = {'tallies':False}
settings.seed = 1
settings.batches = 200
settings.inactive = 50
settings.particles = 1000
settings.source = openmc.source.Source(space=uniform_dist)
settings.temperature = {'multipole': True, 'method': 'interpolation', 'range': [290,2501]}
settings.export_to_xml()

#######################
# MGXS Settings
######################

# ################# Mesh #################
# Create a mesh over the geometry.
# If domain = 'mesh'
# the mgxs will be generated on this mesh.
# If quantify_error = True,
# a 2D fission RR error plot will be generated on this mesh.
mesh = openmc.RegularMesh()
mesh.dimension = [17,17]
mesh.lower_left = [-14,-16.165]
mesh.upper_right = [14,16.165]

############ Set the domain ###########
# If domain_type = 'material',
# assign the list of desired materials to domain (eg [fuel_mat, water_mat]);
# if all materials are desired, use geom.get_all_materials().values()

# If domain_type = 'cell'
# assign the list of desired materials to domain (eg [fuel_cell, water_cell]);
# if all cells are desired, use geom.get_all_material_cells().values()

# If domain_type = 'mesh'
# assign the list of desired meshes to domain (eg [mesh]);

# subdomain_to_plot will be the subdomain whose mgxs are plotted.
# To define subdomain_to_plot
# use the name of a material or cell
# this requires setting a name when defining this material or cell,
# and if domain_type = 'mesh', enter a blank string '' since that cannot be plotted.

domain_type = 'material'
domain = geom.get_all_materials().values()
subdomain_to_plot = fuel_mat

######################
#end example

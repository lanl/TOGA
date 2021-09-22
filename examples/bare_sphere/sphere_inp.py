import openmc
import matplotlib.pyplot as plt
import openmc.checkvalue as cv

#Make OMC input file

num_shells = int(20)/int(10)
sphere_size = 6.38495

u_root = openmc.Universe(universe_id=17)

for shell_ct in range(1,int(num_shells)+1):
    # Create plutonium metal material
    pu = openmc.Material(name='fuel'+str(shell_ct))
    pu.set_density('sum')
    pu.add_nuclide('Pu239', 3.7047e-02)
    pu.add_nuclide('Pu240', 1.7512e-03)
    pu.add_nuclide('Pu241', 1.1674e-04)
    pu.add_element('Ga', 1.3752e-03)
    if shell_ct == 1: materials = openmc.Materials([pu])
    else: materials.append(pu)
    
materials.export_to_xml()

shell_thickness = float(sphere_size)/float(num_shells)
thick = shell_thickness
sphere = []
for shell_ct in range(1,int(num_shells)+1):
    if shell_ct == int(num_shells): 
        sphere.append(openmc.Sphere(r=float(thick), surface_id = int(shell_ct), boundary_type='vacuum'))
    elif shell_ct == 1:
        sphere.append(openmc.Sphere(r=float(thick), surface_id = int(shell_ct)))
    else:
        sphere.append(openmc.Sphere(r=float(thick), surface_id = int(shell_ct)))
    thick += shell_thickness

cells = []
for shell_ct in range(1,int(num_shells)+1):
    if shell_ct == 1:
        cells.append(openmc.Cell(fill=materials[shell_ct-1], region=-sphere[shell_ct-1], cell_id=shell_ct))
    else:
        cells.append(openmc.Cell(fill=materials[shell_ct-1], region=-sphere[shell_ct-1] & +sphere[shell_ct-2], cell_id=shell_ct))
    u_root.add_cells([cells[shell_ct-1]])
geom = openmc.Geometry(u_root)
geom.export_to_xml()

# Finally, define some run settings
settings = openmc.Settings()
settings.batches = 100
settings.inactive = 10
settings.particles = int(1e6)
settings.confidence_intervals = False
settings.export_to_xml()

mesh = openmc.RegularMesh()
mesh.dimension = [170,170]
mesh.lower_left = [-14,-16.165]
mesh.upper_right = [14,16.165]

domain_type = 'material'
domain = geom.get_all_materials().values()
subdomain_to_plot = pu


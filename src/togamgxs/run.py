

import openmc

mpicommand = 'mpirun'
numprocs = '16'
machinefile = 'mach'
numthreads = '7'

openmc.run(mpi_args=[mpicommand,'--bind-to','socket','-np',numprocs,'-machinefile',machinefile],threads=numthreads,cwd='.')

from mpi4py import MPI
from lammps import lammps
import sys

"""
Usage: 
export OMP_NUM_THREADS=<n_openmp_threads>
mpirun -np <n_mpi_processes> python3 run-parallel.py <LAMMPS input file>
"""

input_file = sys.argv[1]

lmp = lammps()
lmp.file(input_file)

me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

print("Proc %d out of %d procs has" % (me,nprocs),lmp)

MPI.Finalize()
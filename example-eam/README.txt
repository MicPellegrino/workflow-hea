I haven't found a way to run in parallel via Jupyter notebook.

From the terminal:

```
export OMP_NUM_THREADS=<n_openmp_threads>
mpirun -np <n_mpi_processes> python3 run-parallel.py <.in file>
```

P.S. the proper way to run with OpenMP is via the OPENMP 
LAMMPS package. However, I am pretty sure it won't work with
Python, since it's not thread safe.
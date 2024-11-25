# workflow-hea

Workflow for simulating HEAs using LAMMPS Python API.

## Initialization

It is suggested to create a virtual environment before starting to run and analyze MD simulations. It is *highly* suggested to use `mamba` for this purpose.
The file `environment.yml` contains an example of a list of modules that you may require. Feel free to modify it, at your own risk.

Create a virtual environment by running:
```
mamba env create --name workflow-hea --file environment.yml
```
activate it by running:
```
mamba activate workflow-hea
```
and deactivate it by running:
```
mamba deactivate
```

I haven't found a way to install Ovito from the `.yml` file that doesn't break MACE. If you want to use Ovito, run the following line after activating the environment:
```
mamba install --strict-channel-priority -c https://conda.ovito.org -c conda-forge ovito=3.11.1
```

## Example: water-graphene simulation

The simulation is described in `water-graphene/SIMULATION-README.txt`. You can use the ´stub´ notebook `water-graphene/workflow.ipynb` to initialize and run a very short simulation, as the full simulation requires either substantial compute resources or a lot of time.

## Example: embedded atom model

The simulation is essentially based on an example by S. Gravelle: https://github.com/simongravelle/lammps-input-files/tree/main/inputs/melting-gold. You can use the notebook `workflow-hea/workflow.ipynb` to initialize and simulate a cube of melting gold using the EAM potential.

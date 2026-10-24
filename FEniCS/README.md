# FEniCS simulations and defect tracking

This directory contains the original finite-element simulation and post-processing code used to study defect-mediated relaxation of a two-dimensional nematic around circular (annular) and elliptical inclusions.

**The simulations have not been rerun in preparing this repository.** The archived outputs used to generate manuscript Figures 13–17 are supplied as processed CSV files in [`../PythonFigs/data/`](../PythonFigs/data/); those figures can be regenerated without FEniCS using [`../PythonFigs/Figures13_17.ipynb`](../PythonFigs/Figures13_17.ipynb).

## Files

| File | Role |
| --- | --- |
| `Annulus/sim.py` | Evolves the relaxational Landau–de Gennes Q-tensor on a circular annulus and writes time-dependent XDMF solutions. |
| `Annulus/data_gen.py` | Reads annular XDMF time series, estimates defect coordinates from a computed diagnostic, and exports a CSV of tracked trajectories. |
| `Ellipse/sim.py` | Corresponding Q-tensor solver for a region between confocal ellipses. |
| `Ellipse/data_gen.py` | Extracts defect trajectories and diagnostic quantities from elliptical XDMF time series and writes a CSV. |

The post-processing scripts also calculate auxiliary quantities from the saved fields, including a free-energy diagnostic.

## Software

The finite-element computations for the manuscript used **FEniCS 2019.1.0** (the legacy FEniCS/DOLFIN interface), together with **mshr** for mesh generation. This source is written for that interface, **not FEniCSx/DOLFINx**.

The simulation scripts also import NumPy and Matplotlib. The post-processing scripts use **meshio**, **PyVista**, **pandas**, and **NumPy**. These simulation and processing dependencies are distinct from the lighter `PythonFigs/requirements.txt`, which is sufficient for plotting the archived CSV data.

An independently tested installation recipe and pinned package environment for rerunning these simulations are not provided here.

## Workflow and paths

1. The `sim.py` files solve the time-dependent Q-tensor equations on the relevant domain. They use spatial finite elements, time stepping, and prescribed boundary orientation, as described in the manuscript.
2. Each solver writes `.xdmf` time-series output (with associated data) to a relative `solutions/` directory. The directory must exist before a run; these large simulation outputs are **not included** in this repository.
3. The corresponding `data_gen.py` reads files from its relative `solutions/` directory, identifies positive and negative defect candidates through extrema of the tracking diagnostic, and records their positions and other quantities in CSV form.
4. The cleaned `PythonFigs/Figures13_17.ipynb` reads the archived CSV data and creates the manuscript plots. The Figure 13 classification also uses a separate processed critical-separation dataset and a Mathematica theoretical curve, both supplied under `PythonFigs/data/`.

The scripts assume relative paths, so working from within the relevant `Annulus/` or `Ellipse/` directory is important if rerunning them.

## Run parameters and archived data

The author adjusted simulation parameters manually between runs. As a result, **the parameters currently hard-coded in the solver files are not a complete record of the settings for each published figure**. In particular:

- `Annulus/sim.py` currently contains an initial-angle choice different from some of the circular-trajectory datasets.
- `Ellipse/sim.py` currently sets `u1 = 1.5`, whereas the Figure 17 processed dataset uses **`u1 = 0.3`** and **`u2 = 2`**.
- The solver scripts select an inclusion-rotation case by a command-line index (`sys.argv[1]`) into an `alpha0_set` array; the arrays themselves can also be changed between runs.
- `Annulus/data_gen.py` currently writes `fig7_defects_data.csv`, whereas the archived plotting dataset is named `calculated_defects_data.csv`. `Ellipse/data_gen.py` writes `ellipse_defects_data.csv`. The processing scripts' current selection arrays and output filenames are not necessarily those used for every archived dataset.

These discrepancies are consistent with the manual run-by-run parameter changes; they should not be interpreted as differing definitions of the governing model. To reproduce any particular raw simulation, first set the desired geometry, initial defect positions, rotation angle, mesh/time settings, and file paths as appropriate for that run.

## Scope of reproducibility

- **Verified separately:** Manuscript Figures 13–17 can be regenerated from the supplied CSV data using the Jupyter plotting notebook.
- **Provided as research source:** FEniCS solvers and post-processing scripts that produced and processed the underlying trajectories.
- **Not supplied or retested:** The full time-dependent XDMF/HDF5 simulation outputs, a complete run-by-run parameter manifest, and a fresh execution of the FEniCS solvers or data-extraction pipeline.

The data filtering and final styling applied to the manuscript figures are documented in `PythonFigs/Figures13_17.ipynb`; they should not be confused with the defect-identification calculations in `data_gen.py`.

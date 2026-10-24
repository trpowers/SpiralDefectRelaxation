# Defect-mediated relaxation of spiral states in frustrated two-dimensional nematics

This repository contains the analytical figure scripts, finite-element simulation source code, and processed numerical data used in the manuscript *Defect-mediated relaxation of spiral states in frustrated two-dimensional nematics*.

The work studies relaxation of a spiral director texture through transport and annihilation of a pair of nematic defects, for circular and elliptical inclusions. The analytical calculations use complex analysis and free-energy minimization; numerical trajectories are obtained from relaxational Landau–de Gennes simulations.

## Repository layout

```text
MathematicaFigs/
    README.md
    Fig1.wl, Fig2.wl, ..., Fig6.wl
    Figs7_8_9.wl
    Fig10.wl, Fig11.wl, Fig12.wl

PythonFigs/
    README.md
    Figures13_17.ipynb
    requirements.txt
    data/
        calculated_defects_data.csv
        critical_phi_calculated_defects_data.csv
        ellipse_defects_data.csv
        fig7data_mathematica.csv
        fig9_mathematica_data.csv

FEniCS/
    README.md
    Annulus/
        sim.py
        data_gen.py
    Ellipse/
        sim.py
        data_gen.py
```

- **[MathematicaFigs](MathematicaFigs/README.md):** Wolfram Language scripts that generate Figures 1–12 as PDFs. Figures 7–9 share one script. The manuscript calculations used Mathematica 13.3.1; selected scripts also require MaTeX and a LaTeX installation.
- **[PythonFigs](PythonFigs/README.md):** A Jupyter notebook that regenerates Figures 13–17 directly from the five supplied CSV files. It does **not** rerun FEniCS or Mathematica. The notebook exports both PDF and PNG figures.
- **[FEniCS](FEniCS/README.md):** Original FEniCS solver and post-processing scripts for annular and elliptical domains. These document the simulation methodology but have not been rerun or independently validated as part of the repository preparation.

## Reproducing the figures

**Figures 1–12:** Follow [`MathematicaFigs/README.md`](MathematicaFigs/README.md) to execute the corresponding `.wl` scripts in Mathematica.

**Figures 13–17:** From `PythonFigs/`, install the notebook dependencies and open the notebook:

```bash
python -m pip install -r requirements.txt
jupyter notebook Figures13_17.ipynb
```

Run the notebook cells in order. The required CSV files are included under `PythonFigs/data/`. Output PDFs and PNGs are written to `PythonFigs/figures/`. **FEniCS is not required to regenerate Figures 13–17 from the archived data.**

## Simulation-code provenance

The FEniCS scripts were used during the research with parameter values changed manually between simulation runs. Their current hard-coded defaults therefore do not necessarily match the settings used for every archived dataset. The original time-dependent solution files are not included here; the processed CSV data used for the published figures are provided in `PythonFigs/data/`.

For the FEniCS model, dependencies, post-processing workflow, and limitations, see [`FEniCS/README.md`](FEniCS/README.md).

## Manuscript

*Defect-mediated relaxation of spiral states in frustrated two-dimensional nematics* (submitted to *Soft Matter*). Publication information can be added here when available.

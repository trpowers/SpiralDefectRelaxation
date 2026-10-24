# Figures 13–17 (processed finite-element data)

`Figures13_17.ipynb` regenerates manuscript Figures 13–17 from the five CSV files in `data/`. No FEniCS simulations or Mathematica calculations are run. Each figure is in a separate notebook cell for editing.

## Running

```bash
python -m pip install -r requirements.txt
jupyter notebook Figures13_17.ipynb
```

Run all cells in order. Each figure is displayed in the notebook and saved as a vector PDF plus a 300-dpi PNG in the automatically created `figures/` directory. `data/` must be alongside the notebook. Nothing depends on machine-specific paths.

The notebook is saved without cached outputs to keep the repository compact and easy to maintain. Run the cells to regenerate the figures. `figures/` is gitignored by default; remove that line from `.gitignore` if you want to commit generated PDFs.

## Reproduction notes

- Figure 13: `max_delta >= 50`, `min_delta <= -50`, 55 phase markers, and the supplied 180-point Mathematica CSV.
- Figure 14: annulus tracking cutoff 50; alpha = 0.6 pi; initial `phi_param = 0.52` rad.
- Figure 15: annulus tracking cutoff 10; alpha = 0.8 pi.
- Figure 16: annulus tracking cutoff 10; initial `phi_param = 0.26` rad; six alpha series.
- Figure 17: ellipse tracking cutoff 100, `t > 0`, `u1 = 0.3`, `u2 = 2`, and six alpha series. 

The notebook changes implementation (one coherent set of imports, relative paths, compact cells, PDF+PNG export) but does not rerun FEM or independently validate those underlying simulation results.

Fonts: Times New Roman is selected when installed, otherwise STIX Two Text or STIXGeneral. Appearance may vary slightly between operating systems.

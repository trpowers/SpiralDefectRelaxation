# Mathematica figures (Figures 1–12)

These Wolfram Language (`.wl`) scripts generate the analytical and director-field figures for the *Soft Matter* manuscript on defect-mediated relaxation of nematic spiral textures. The figures are exported as PDF files. The finite-element simulation figures (13–17) are generated separately from saved CSV data using Python/Jupyter.

## Requirements

- **Wolfram Mathematica 13.3.1**, the version used for the manuscript calculations.
- **MaTeX** and a working LaTeX installation for Figures **2, 6, 7–9, and 12**, which use `Needs["MaTeX`"]` for mathematical labels. MaTeX must be installed and configured before running those scripts.
- No external input data files are read by these scripts.

## Files

| Script | PDF output | Contents |
| --- | --- | --- |
| `Fig1.wl` | `Fig1.pdf` | Defect-free spiral director textures in a circular annulus |
| `Fig2.wl` | `Fig2.pdf` | Conformal maps of the annulus to a strip and upper half-plane, with defect director fields |
| `Fig3.wl` | `Fig3.pdf` | Half-integer defects and reflected images |
| `Fig4.wl` | `Fig4.pdf` | Defect-pair transport and change of spiral charge |
| `Fig5.wl` | `Fig5.pdf` | Defect-pair director textures at different inclusion rotations and separations |
| `Fig6.wl` | `Fig6.pdf` | Integration contours surrounding a defect pair |
| `Figs7_8_9.wl` | `Fig7.pdf`, `Fig8.pdf`, `Fig9.pdf` | Circular-annulus optimal radius, free-energy barriers, and critical angular separation |
| `Fig10.wl` | `Fig10.pdf` | Defect-free spiral textures for an elliptic inclusion |
| `Fig11.wl` | `Fig11.pdf` | Defect-mediated unwinding for an elliptic inclusion |
| `Fig12.wl` | `Fig12.pdf` | Free-energy curves and barrier locations for the needle-like inclusion |

There is one shared script for Figures 7–9 because those plots use common numerical calculations.

## Running a script

1. Open Mathematica and start a kernel. For reliable execution, use a fresh kernel when switching between scripts, since several define or clear symbols in the global context.
2. Set the working directory to the folder containing the `.wl` files. For example, in a notebook saved in that folder:

   ```wolfram
   SetDirectory[NotebookDirectory[]];
   ```

   Alternatively, use `SetDirectory["/path/to/Mathematica"]` with the actual path on your machine.
3. Load a script with `Get`, for example:

   ```wolfram
   Get["Fig2.wl"];
   ```

   To generate Figures 7–9, run `Get["Figs7_8_9.wl"];` once.

PDFs are written to **Mathematica's current working directory** (`Directory[]`), not automatically to a separate `figures/` folder. You can change this by editing the `Export` commands near the ends of the scripts.

## Reproducibility notes

- Parameters, figure styling, and numerical settings are defined within the scripts. The energy plots use numerical optimization and interpolation; their execution may take longer than the schematic director plots.
- The scripts generate the PDFs named above. The filenames in the manuscript's `figures/` folder may use older descriptive names (for example, `FtotWithBarriersNeedle.pdf` instead of `Fig12.pdf`). Rename or copy the generated PDFs as necessary when updating the manuscript.
- This folder contains the analytical figure scripts only. It does **not** rerun the FEniCS simulations underlying Figures 13–17.

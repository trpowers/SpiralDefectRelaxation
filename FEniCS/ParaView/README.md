# ParaView visualization state

This directory contains `pview_state.pvsm`, a ParaView visualization-state file associated with the time-dependent finite-element calculations in **Defect-mediated relaxation of spiral states in frustrated two-dimensional nematics**.

The state was saved in **ParaView 6.0.1** (the version recorded in its XML). It preserves a visualization pipeline, display properties, a camera, and animation settings. It **does not contain the underlying simulation fields or a finished video**.

## What this state describes

The state references **one circular-annulus simulation**, not all eight supplementary videos. Its data source is an XDMF file named:

```text
defects_solution_alpha01.57_phi0.52_r0.61_KA0.0010_rho0.37.xdmf
```

The parameters encoded in that filename are for **Video S2**, a circular example near `alpha = pi/2`, `rho = 1/e`, and with an initial angle near `pi/6`. 

The saved pipeline includes:

- An **Xdmf3ReaderS** reader for the time-dependent XDMF solution.
- Point-data arrays named `Qxx`, `Qxy`, `S_order_param`, and `director`.
- A glyph filter stored as `an_glyph`, oriented by the `director` array, with saved scale factor `0.05` and a uniform spatial sampling mode.
- A surface representation colored by `S_order_param` and a time annotation displaying the current time in units of `gamma/A`.
- Camera/view settings and an animation scene. The XDMF reader lists 1,201 time values from 0 to 1199.01; the state was last saved at time 1003.01.

These values describe the **saved state**, not necessarily the export settings or complete time interval used for the published movie.

## How to open the file

1. Open **ParaView**. Version **6.0.1** is the closest match to the state; other versions may require adjustments.
2. Choose **File > Load State...** and select `pview_state.pvsm`. This is *not* opened using the ordinary **File > Open** command.
3. In the **Load State Options** dialog, choose **Choose File Names** and select the matching `.xdmf` solution file on your machine. Alternatively, select **Search files under specified directory** and browse to the directory containing it.
4. Make sure the `.xdmf` file's associated field-data files (if externally referenced by XDMF) are present in their expected locations.
5. Inspect the restored pipeline and use the animation controls to advance through its time steps. To make a movie, use ParaView's animation-export function and select the desired frame range, frame rate, resolution, and output format.

The original state contains an **absolute path on the author's computer**. That path will not work on another machine; select the replacement data file when loading the state.

### Possible custom-filter dependency

The file stores a filter with internal type **`an_glyph`**. It is not identified within the state as a built-in filter or as a particular external plugin. If ParaView reports an unknown proxy/filter when loading, the visualization may require a custom filter or plugin from the original environment. Ask the simulation author for its definition or an equivalent standard-ParaView pipeline. The state file does not provide that implementation.

## Reproducibility and scope

- The `.pvsm` file is **visualization metadata**, not a numerical solution. The processed defect-position CSVs in `PythonFigs/data/` cannot replace the full spatial XDMF solution needed to recreate the animation.
- The associated XDMF dataset and any companion arrays are **not included here**; nor are the finished AVI files. The supplementary movies are supplied separately with the journal article.
- This one state does **not establish the visualization settings for every Video S1–S8**. Other runs may require changing the source file, its geometry, and visualization settings, or loading additional state files.
- Video encoding, actual AVI frame rate, final output resolution, and exact exported time interval **cannot be determined reliably from this state**.
- This README was prepared by inspecting the state XML. The visualization has **not been independently restored or used to export a movie**, because the referenced full-field simulation dataset was not supplied.

For background on state-file loading, see the [ParaView documentation on saving and loading states](https://docs.paraview.org/en/latest/UsersGuide/savingResults.html).

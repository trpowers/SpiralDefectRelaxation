(* ::Package:: *)

(* ============================================================ *)
(* FIGURE 1: Defect-free spiral textures                        *)
(* ============================================================ *)

(* Run this script from Mathematica using:
     SetDirectory["/path/to/working/directory"];
     Get["Fig1.wl"];

   Alternatively, from a notebook in the same directory:
     SetDirectory[NotebookDirectory[]];
     Get["Fig1.wl"];

   The script exports Fig1.pdf to the current working directory,
   as returned by Directory[].
   Requires no external packages or input files.
*)



(* Figure 1: Defect-free spiral director fields in a circular annulus.

   Geometry: outer radius 1, inner radius rho = Exp[-L].
   Boundary conditions: theta(1) = 0, theta(rho) = alpha.
   Spiral solution: theta(r) = alpha Log[r]/Log[rho].
   With an outward radial integration path, spiral charge mu = -alpha.

   The four panels use alpha = 0, Pi/2, Pi, 3 Pi/2.
*)

ClearAll["Global`*"];

(* Parameters *)
L = 1;
rho = Exp[-L];
alphaVals = {0, Pi/2, Pi, 3 Pi/2};

(* Draw streamlines of the unit director for a given inner boundary angle. *)
makePlot[alpha0_] := Module[{thetaSp, director, fieldPlot},
  thetaSp = alpha0 Log[Sqrt[x^2 + y^2]]/Log[rho];
  director = {Cos[thetaSp], Sin[thetaSp]};

  fieldPlot = StreamPlot[
    director, {x, -1, 1}, {y, -1, 1},
    RegionFunction -> Function[{x, y, vx, vy},
      rho <= Sqrt[x^2 + y^2] <= 1],
    StreamMarkers -> None,
    StreamColorFunction -> None,
    RegionBoundaryStyle -> None,
    RegionFillingStyle -> None,
    Frame -> False, Axes -> False, AspectRatio -> 1,
    PlotRange -> {{-1, 1}, {-1, 1}},
    PlotRangePadding -> None,
    ImagePadding -> 0, ImageSize -> 250
  ];

  Show[
    fieldPlot,
    Graphics[{
      White, Disk[{0, 0}, rho + 0.002],
      Black, Thick, Circle[{0, 0}, rho], Circle[{0, 0}, 1]
    }],
    PlotRange -> 1.04 {{-1, 1}, {-1, 1}},
    Frame -> False, Axes -> False, AspectRatio -> 1,
    ImagePadding -> 5
  ]
];

(* Four-panel figure *)
plt1 = GraphicsRow[makePlot /@ alphaVals, Spacings -> 20];

(* Save a vector PDF in the current working directory. *)
outputFile = FileNameJoin[{Directory[], "Fig1.pdf"}];
Export[outputFile, plt1, "PDF", ImageResolution -> 150];
Print["Figure 1 saved to: ", outputFile];

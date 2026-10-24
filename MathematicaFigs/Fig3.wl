(* ============================================================ *)
(* FIGURE 3: Half-integer defects and their reflected images     *)
(* ============================================================ *)
(* Produces Fig3.pdf in Mathematica's current working directory. *)

ClearAll["Global`*"];

(* Geometry: the real defect lies in the upper half-plane; its image
   is reflected across the real axis y = 0. *)
realDefect = {0, 1};
imageDefect = {0, -1};
H = 3;

(* Director angles for a defect together with its image.

   The real and image +1/2 defects have orientation constants
   +Pi/4 and -Pi/4. For the -1/2 pair these are +3 Pi/4 and
   -3 Pi/4. In both cases the constants cancel in the total angle.

   Mathematica's two-argument ArcTan[dx, dy] returns the polar angle
   Arg[dx + I dy]. Thus thetaPlus is half the sum of the polar angles
   measured from the real defect and its image; thetaMinus is its
   negative. On y = 0 these two angles cancel, so theta = 0 there. *)

thetaPlus[x_, y_] :=
  1/2 (ArcTan[x - realDefect[[1]], y - realDefect[[2]]] +
       ArcTan[x - imageDefect[[1]], y - imageDefect[[2]]]);

thetaMinus[x_, y_] := -thetaPlus[x, y];

nfieldPlus[x_, y_] := {Cos[thetaPlus[x, y]], Sin[thetaPlus[x, y]]};
nfieldMinus[x_, y_] := {Cos[thetaMinus[x, y]], Sin[thetaMinus[x, y]]};

(* Use the original StreamPlot settings and distinguish real defects
   (opaque) from their images (partially transparent). *)
makePlot[nfield_, markerColor_] :=
  StreamPlot[
    nfield[x, y], {x, -H, H}, {y, -H, H},
    StreamMarkers -> None,
    StreamColorFunction -> None,
    Frame -> False,
    Axes -> False,
    AspectRatio -> 1,
    PlotRange -> {{-H, H}, {-H, H}},
    Epilog -> {
      markerColor, Disk[realDefect, 0.08],
      Opacity[0.35], markerColor, Disk[imageDefect, 0.08],
      Opacity[1], Black, AbsoluteThickness[3],
      Line[{{-H, 0}, {H, 0}}]
    },
    ImageSize -> 350
  ];

fig3 = GraphicsRow[
  {makePlot[nfieldPlus, Blue], makePlot[nfieldMinus, Red]},
  Spacings -> 25
];

Export["Fig3.pdf", fig3, "PDF", ImageResolution -> 150];

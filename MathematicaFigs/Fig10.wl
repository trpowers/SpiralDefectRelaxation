(* Figure 10: Defect-free spiral states in a confocal elliptic annulus.

   Coordinates: z = Exp[I alpha] Cosh[u + I v], with u1 <= u <= u2.
   The focal half-distance is the unit of length. Both boundary ellipses
   rotate through alpha, but the outer director remains along the fixed
   laboratory x-axis. This is the geometry used in the manuscript.

   Model: one-constant Frank elasticity, with theta = alpha on u = u1
   and theta = 0 on u = u2. The parameters must satisfy 0 < u1 < u2.

   Dependencies: Mathematica only; MaTeX is not required.
   Output: Fig10.pdf in Mathematica's current working directory,
   Directory[], which need not be the directory containing this script.
*)

ClearAll[x, y, t, u0, \[Alpha], u1, u2, \[Alpha]Vals, u, theta,
  ellipsePoint, makePlot, plt1];

(* Geometry and accumulated counterclockwise rotations.
   The conformal strip width is ell = u2 - u1 = 1.9. *)
u1 = 0.1;
u2 = 2;
\[Alpha]Vals = {0, Pi/2, Pi, 3 Pi/2};

(* Inverse elliptic radial coordinate in the rotated geometry.
   Although the principal ArcCosh has a branch cut, its real part is
   smooth throughout this annulus. The focal segment lies inside the
   excluded inclusion; the remaining cut changes only the imaginary part. *)
u[x_, y_, \[Alpha]_] := Re[ArcCosh[Exp[-I \[Alpha]] (x + I y)]];

(* Director angle measured in the laboratory frame.
   This is linear in u and therefore harmonic in the one-constant model.
   Keep the full alpha, not alpha modulo Pi: it distinguishes the spiral
   sectors. The inner-to-outer spiral charge is mu = -alpha. *)
theta[x_, y_, \[Alpha]_] :=
  \[Alpha] (u2 - u[x, y, \[Alpha]])/(u2 - u1);

(* Cartesian form of Exp[I alpha] Cosh[u0 + I t].
   An ellipse at u0 has semiaxes Cosh[u0] and Sinh[u0], and the
   laboratory-frame foci are at z = +/- Exp[I alpha]. *)
ellipsePoint[u0_, t_, \[Alpha]_] := {
  Cosh[u0] Cos[t] Cos[\[Alpha]] - Sinh[u0] Sin[t] Sin[\[Alpha]],
  Cosh[u0] Cos[t] Sin[\[Alpha]] + Sinh[u0] Sin[t] Cos[\[Alpha]]
};

(* One panel. The original StreamPlot, boundary, and layout settings
   are preserved, including the 8% margin around the outer semimajor axis. *)
makePlot[\[Alpha]_] :=
  Module[{nfield, streamPlot, boundaryPlot, rPlot},
    nfield = {Cos[theta[x, y, \[Alpha]]], Sin[theta[x, y, \[Alpha]]]};
    rPlot = 1.08 Cosh[u2];

    streamPlot = StreamPlot[
      nfield, {x, -rPlot, rPlot}, {y, -rPlot, rPlot},
      RegionFunction ->
        Function[{x, y, vx, vy}, u1 <= u[x, y, \[Alpha]] <= u2],
      StreamMarkers -> None, StreamColorFunction -> None,
      RegionBoundaryStyle -> None, RegionFillingStyle -> None,
      Frame -> False, Axes -> False, AspectRatio -> 1,
      PlotRange -> {{-rPlot, rPlot}, {-rPlot, rPlot}},
      PlotRangePadding -> None, ImagePadding -> 0, ImageSize -> 250
    ];

    boundaryPlot = ParametricPlot[
      Evaluate[ellipsePoint[#, t, \[Alpha]] & /@ {u1, u2}],
      {t, 0, 2 Pi}, PlotStyle -> Directive[Black, Thick],
      Axes -> False, Frame -> False
    ];

    Show[streamPlot, boundaryPlot,
      PlotRange -> {{-rPlot, rPlot}, {-rPlot, rPlot}},
      Frame -> False, Axes -> False, AspectRatio -> 1, ImagePadding -> 5
    ]
  ];

(* Four panels, in the same order and with the original spacing. *)
plt1 = GraphicsRow[makePlot /@ \[Alpha]Vals, Spacings -> 20];

Export["Fig10.pdf", plt1, "PDF", ImageResolution -> 150];

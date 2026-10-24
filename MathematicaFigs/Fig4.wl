(* ================================================================ *)
(* FIGURE 4: Defect-pair transport and the change of spiral charge   *)
(* ================================================================ *)

(* Reproduces the two-row, four-column figure with the real and
   image-induced defect director fields, a sector contour, and the
   branch cut of the principal complex logarithm.

   Geometry: rho < |z| < 1, L = Log[1/rho].
   Defects: zPlus = rPlus Exp[I phi], zMinus = rMinus Exp[-I phi].
   Spiral angle: thetaSp = alpha Log[|z|]/Log[rho].

   Uses Mathematica's built-in QPochhammer (no external packages).
   Exports Fig4.pdf to the current working directory (Directory[]).
*)

ClearAll[
  rho, a, b, H, L, q, rPlus, rMinus, Zmap, Pcrowdy,
  phi1Contour, phi2Contour, deltaContour, rInnerContour,
  rOuterContour, contourPts, makeDefectPlot, phiVals, pltGrid
];

(* Annulus and prime-function parameter. *)
rho = 1/E;
a = rho;
b = 1;
H = b;
L = Log[1/rho];
q = Exp[-Pi^2/L];

rPlus = 0.7;
rMinus = 0.7;

(* Manuscript maps: Z = Log[z/rho], zeta = -Exp[-I Pi Z/L].
   Since rho > 0 and Log[rho] = -L, we have Z = Log[z] + L, so
   -Exp[-I Pi (Log[z] + L)/L] = Exp[-I Pi Log[z]/L].
   The two minus signs cancel because Exp[-I Pi] = -1.
   Principal Log[z] retains the original branch-cut convention. *)
Zmap[z_] := Exp[-I Pi Log[z]/L];

(* Exact Schottky-Klein prime-function product (no truncation). *)
Pcrowdy[lambda_] :=
  QPochhammer[lambda, q^2] QPochhammer[q^2/lambda, q^2];

(* Closed magenta sector contour. Its sampling is deliberately fixed. *)
phi1Contour = 0.9 Pi/3;
phi2Contour = Pi;
deltaContour = 0.02;
rInnerContour = a + deltaContour;
rOuterContour = b - deltaContour;

contourPts = Join[
  Table[
    rOuterContour {Cos[t], Sin[t]},
    {t, phi1Contour, phi2Contour, (phi2Contour - phi1Contour)/60}
  ],
  Table[
    s {Cos[phi2Contour], Sin[phi2Contour]},
    {s, rOuterContour, rInnerContour, -(rOuterContour - rInnerContour)/30}
  ],
  Table[
    rInnerContour {Cos[t], Sin[t]},
    {t, phi2Contour, phi1Contour, -(phi2Contour - phi1Contour)/60}
  ],
  Table[
    s {Cos[phi1Contour], Sin[phi1Contour]},
    {s, rInnerContour, rOuterContour, (rOuterContour - rInnerContour)/30}
  ]
];

(* One annulus panel, for defect half-separation phi and spiral angle alpha. *)
makeDefectPlot[phi_, alpha_ : 0] := Module[
  {
    zPlus, zMinus, zetaPlus, zetaMinus,
    Wcut, twoIOmega, Omega, rawTheta, thetaPlot, nfieldDef,
    branchCutContour, branchCutLines, stream, decorations
  },

  zPlus = rPlus Exp[I phi];
  zMinus = rMinus Exp[-I phi];
  zetaPlus = Zmap[zPlus];
  zetaMinus = Zmap[zMinus];

  (* Prime-function quotient W. It is the argument of the logarithm
     in the defect complex potential. Keeping W itself accessible
     also locates the principal-Log branch cut. *)
  Wcut[x_?NumericQ, y_?NumericQ] := Module[
    {z = x + I y, zeta, lambdaPlus, lambdaMinus,
     lambdaTildePlus, lambdaTildeMinus},
    zeta = Zmap[z];
    lambdaPlus = zetaPlus/zeta;
    lambdaMinus = zetaMinus/zeta;
    lambdaTildePlus = Conjugate[zetaPlus]/zeta;
    lambdaTildeMinus = Conjugate[zetaMinus]/zeta;
    (Pcrowdy[lambdaPlus] Pcrowdy[lambdaTildePlus])/
      (Pcrowdy[lambdaMinus] Pcrowdy[lambdaTildeMinus])
  ];

  (* On the principal logarithm's branch cut, Im[W] = 0 and Re[W] < 0.
     This numerical contour is deliberately preserved, including
     its original resolution and the annular region restriction. *)
  branchCutContour = ContourPlot[
    Im[Wcut[x, y]], {x, -H, H}, {y, -H, H},
    Contours -> {0},
    ContourShading -> False,
    RegionFunction -> Function[{x, y, z},
      a <= Sqrt[x^2 + y^2] <= b && Re[Wcut[x, y]] < 0
    ],
    ContourStyle -> Black,
    PlotPoints -> 100,
    MaxRecursion -> 3,
    Frame -> False,
    Axes -> False,
    BoundaryStyle -> None
  ];
  branchCutLines = Cases[
    Normal[branchCutContour], Line[p_] :> Line[p], Infinity
  ];

  (* Defect angle: Re[-I Log[W]/2] = Arg[W]/2.
     Retain the complex-potential convention of the original code. *)
  twoIOmega[z_?NumericQ] := Log[Wcut[Re[z], Im[z]]];
  Omega[z_?NumericQ] := -I/2 twoIOmega[z];
  rawTheta[x_?NumericQ, y_?NumericQ] := Re[Omega[x + I y]];
  thetaPlot[x_?NumericQ, y_?NumericQ] :=
    rawTheta[x, y] + alpha Log[Sqrt[x^2 + y^2]]/Log[rho];
  nfieldDef[x_?NumericQ, y_?NumericQ] :=
    {Cos[thetaPlot[x, y]], Sin[thetaPlot[x, y]]};

  (* Original streamline options are retained. *)
  stream = StreamPlot[
    Evaluate[nfieldDef[x, y]], {x, -H, H}, {y, -H, H},
    RegionFunction -> Function[{x, y, vx, vy},
      a <= Sqrt[x^2 + y^2] <= b
    ],
    RegionFillingStyle -> None,
    RegionBoundaryStyle -> None,
    StreamMarkers -> None,
    StreamColorFunction -> None,
    StreamPoints -> Medium,
    StreamScale -> {Tiny, 1.2, None},
    StreamStyle -> Directive[AbsoluteThickness[1.2]],
    PerformanceGoal -> "Quality",
    Frame -> False,
    Axes -> False,
    AspectRatio -> 1,
    PlotRange -> 1.02 {{-H, H}, {-H, H}},
    PlotRangePadding -> None,
    ImagePadding -> 5,
    ImageSize -> 220
  ];

  (* Sector contour, arrow, principal-Log branch cut, and defects. *)
  decorations = Graphics[{
    White, Disk[{0, 0}, a + 0.002],
    Black, Thick, Circle[{0, 0}, a], Circle[{0, 0}, b],
    Magenta, Thick, Line[contourPts],
    Arrowheads[0.1],
    Arrow[{0.85 {Cos[Pi], Sin[Pi]}, 0.45 {Cos[Pi], Sin[Pi]}}],
    GrayLevel[0.35], AbsoluteThickness[3], Dashing[{0.04, 0.025}],
    branchCutLines,
    Blue, Disk[{Re[zPlus], Im[zPlus]}, 0.05],
    Red, Disk[{Re[zMinus], Im[zMinus]}, 0.05]
  }];

  Show[
    stream, decorations,
    PlotRange -> 1.02 {{-H, H}, {-H, H}},
    AspectRatio -> 1,
    ImagePadding -> 5,
    ImageSize -> 220
  ]
];

(* Columns: defect progression. Rows: alpha = 0 and alpha = Pi. *)
phiVals = {Pi/12, Pi/3, 2 Pi/3, 11 Pi/12};
pltGrid = GraphicsGrid[
  {makeDefectPlot[#, 0] & /@ phiVals,
   makeDefectPlot[#, Pi] & /@ phiVals},
  Spacings -> {10, 10}
];

Export["Fig4.pdf", pltGrid, "PDF", ImageResolution -> 150];

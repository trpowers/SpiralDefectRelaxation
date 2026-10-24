(* ============================================================ *)
(* Figure 2: Conformal maps and defect-pair director fields       *)
(* ============================================================ *)
(*
   PURPOSE
   Reproduce the three-panel figure showing the annulus, the
   logarithmic strip, and the upper-half zeta plane, with the
   director field, +/-1/2 defects, periodic copies, and images.

   MANUSCRIPT CONVENTIONS
   rho < |z| < 1, ell = Log[1/rho], z+/- = R Exp[+/- I dphi/2].
   Z = Log[z/rho] = u + I v, where 0 < u < ell.
   zeta = -Exp[-I Pi Z/ell].

   The strip is shown over several periods in v. The rods in
   each panel depict the director angle theta as a scalar
   function of the indicated coordinates; they are not tangent
   vectors transformed by the derivative of the conformal map.

   DEPENDENCIES
   Mathematica, the external MaTeX package, and a LaTeX installation
   with newtxtext/newtxmath packages. No external data files.

   OUTPUT
   Exports Fig2.pdf to Directory[].
*)

ClearAll["Global`*"];

Needs["MaTeX`"];

(* ============================================================ *)
(* Parameters                                                   *)
(* ============================================================ *)

rho = 1/E;
ell = Log[1/rho];

(* Annular nome; machine precision is used for rod sampling. *)
q = Exp[-Pi^2/ell];
qNum = N[q];

alpha = Pi/2;
R = 0.75;
dphi = Pi/4;

(* ============================================================ *)
(* Styling                                                      *)
(* ============================================================ *)

directorColor = RGBColor[0.29, 0.46, 0.6];
directorThickness = AbsoluteThickness[1.0];

plusColor = RGBColor[0.20, 0.35, 0.80];   (* blue: +1/2 *)
minusColor = RGBColor[0.80, 0.20, 0.20];  (* red:  -1/2 *)
clusterColor = RGBColor[0.55, 0.30, 0.75];

innerBoundaryColor = RGBColor[0.10, 0.60, 0.20];  (* green *)
outerBoundaryColor = RGBColor[0.95, 0.55, 0.10];  (* orange *)

guideColor = GrayLevel[0.50];
guideThickness = AbsoluteThickness[1.8];
guideDash = AbsoluteDashing[{4, 3}];

(* ============================================================ *)
(* MaTeX settings (must precede any panel rendering)           *)
(* ============================================================ *)

maTeXPreamble = {
     "\\usepackage{newtxtext}",
     "\\usepackage{newtxmath}"
   };

mapMagnification = 1.25;
planeMagnification = 1.15;

(* ============================================================ *)
(* Figure sizing                                                *)
(* ============================================================ *)

annulusImageSize = 210;
stripImageWidth = 155;
stripAspectRatio = 2.7;
stripEndPad = Pi/3;

uhpImageWidth = 190;
uhpXMax = 4.5;
uhpYMax = 6.0;
uhpBelowAxisFactor = 0.70;

gapAnnulusStrip = 25;
gapStripUHP = 30;

(* ============================================================ *)
(* Rod lengths and sampling density                             *)
(* ============================================================ *)

directorGridScale = 0.5;

annulusRodLength = 0.045;
stripRodLength = 0.04;
uhpRodLength = 0.18;

annulusGridSpacing = 0.18 directorGridScale;
stripUSpacing = (ell/5) directorGridScale;
stripVSpacing = (Pi/6) directorGridScale;
uhpXSpacing = 0.70 directorGridScale;
uhpYSpacing = 0.60 directorGridScale;

(* ============================================================ *)
(* Defect markers                                               *)
(* ============================================================ *)

defectMarkerSize = 5;
clusterMarkerSize = 6;
imageOpacity = 0.75;

(* All point markers are circular, with size specified in printer points.
   Defect markers and the branch-point cluster use the same primitive. *)
pointMarker[pt_, color_, size_, opacity_ : 1] :=
  Inset[
    Graphics[
      {Opacity[opacity], color, Disk[{0, 0}, 1]},
      PlotRange -> {{-1.1, 1.1}, {-1.1, 1.1}},
      ImagePadding -> 0,
      Background -> None,
      ImageSize -> size
    ],
    pt, Center
  ];

defectMarker[pt_, color_, opacity_ : 1] :=
  pointMarker[pt, color, defectMarkerSize, opacity];

clusterMarker[pt_] :=
  pointMarker[pt, clusterColor, clusterMarkerSize];

(* ============================================================ *)
(* Complex-plane corner symbol                                  *)
(* ============================================================ *)

planeMarker[pos_, label_, sx_, sy_, labelFrac_ : {0.52, 0.60}] :=
    {
       Black,
       AbsoluteThickness[1.2],
       Line[{pos + {0, sy}, pos, pos + {sx, 0}}],
       Inset[label, pos + {labelFrac[[1]] sx, labelFrac[[2]] sy}, Center]
     };

(* ============================================================ *)
(* Crowdy prime-function product: nome parameter q^2           *)
(* ============================================================ *)

Pcrowdy[lambda_?NumericQ] :=
    Module[{lam = N[lambda]},
       QPochhammer[lam, qNum^2] * QPochhammer[qNum^2/lam, qNum^2]
     ];

(* ============================================================ *)
(* Defect locations and conformal maps                          *)
(* ============================================================ *)

zPlus[R_?NumericQ, dphi_?NumericQ] := N[R] Exp[I N[dphi]/2];
zMinus[R_?NumericQ, dphi_?NumericQ] := N[R] Exp[-I N[dphi]/2];

(* Real strip coordinate u = Log[|z|/rho], with 0 < u < ell. *)
u0FromR[R_?NumericQ] := Log[N[R]/rho];

zetaFromStrip[u_?NumericQ, v_?NumericQ] :=
    -Exp[-I Pi (N[u] + I N[v])/ell];

zetaPlus[R_?NumericQ, dphi_?NumericQ] :=
    zetaFromStrip[u0FromR[R], N[dphi]/2];

zetaMinus[R_?NumericQ, dphi_?NumericQ] :=
    zetaFromStrip[u0FromR[R], -N[dphi]/2];

(* Image defects are reflections across the real zeta axis. *)
zetaPlusImage[R_?NumericQ, dphi_?NumericQ] :=
  Conjugate[zetaPlus[R, dphi]];
zetaMinusImage[R_?NumericQ, dphi_?NumericQ] :=
  Conjugate[zetaMinus[R, dphi]];

(* ============================================================ *)
(* Director fields                                              *)
(* ============================================================ *)

(* Prime-function quotient for the +/-1/2 defect pair and images.
   The defect director angle is half its complex argument. *)
defectRatioZeta[zeta_?NumericQ, R_?NumericQ, dphi_?NumericQ] :=
    Module[{zp, zm},
       zp = zetaPlus[R, dphi];
       zm = zetaMinus[R, dphi];
       (Pcrowdy[zp/zeta] Pcrowdy[Conjugate[zp]/zeta])/
          (Pcrowdy[zm/zeta] Pcrowdy[Conjugate[zm]/zeta])
     ];

thetaDefZeta[zeta_?NumericQ, R_?NumericQ, dphi_?NumericQ] :=
    1/2 Arg[defectRatioZeta[zeta, R, dphi]];

thetaAnnulus[x_?NumericQ, y_?NumericQ, alpha_?NumericQ, R_?NumericQ,
   dphi_?NumericQ] :=
    Module[{z, r, u, v, zeta},
       z = N[x] + I N[y];
       r = Abs[z];
       u = Log[r/rho];
       v = Arg[z];
       zeta = zetaFromStrip[u, v];
       N[alpha] Log[r]/Log[rho] + thetaDefZeta[zeta, R, dphi]
     ];

nAnnulus[x_?NumericQ, y_?NumericQ, alpha_?NumericQ, R_?NumericQ,
   dphi_?NumericQ] :=
    Module[{theta},
       theta = thetaAnnulus[x, y, alpha, R, dphi];
       {Cos[theta], Sin[theta]}
     ];

(* Log[r]/Log[rho] = 1-u/ell in strip coordinates. *)
thetaStrip[u_?NumericQ, v_?NumericQ, alpha_?NumericQ, R_?NumericQ,
   dphi_?NumericQ] :=
    Module[{zeta},
       zeta = zetaFromStrip[u, v];
       N[alpha] (1 - N[u]/ell) + thetaDefZeta[zeta, R, dphi]
     ];

nStrip[u_?NumericQ, v_?NumericQ, alpha_?NumericQ, R_?NumericQ,
   dphi_?NumericQ] :=
    Module[{theta},
       theta = thetaStrip[u, v, alpha, R, dphi];
       {Cos[theta], Sin[theta]}
     ];

(* In the upper half-plane, 1-u/ell = Arg[zeta]/Pi. *)
thetaUHP[x_?NumericQ, y_?NumericQ, alpha_?NumericQ, R_?NumericQ,
   dphi_?NumericQ] :=
    Module[{zeta},
       zeta = N[x] + I N[y];
       N[alpha] Arg[zeta]/Pi + thetaDefZeta[zeta, R, dphi]
     ];

nUHP[x_?NumericQ, y_?NumericQ, alpha_?NumericQ, R_?NumericQ,
   dphi_?NumericQ] :=
    Module[{theta},
       theta = thetaUHP[x, y, alpha, R, dphi];
       {Cos[theta], Sin[theta]}
     ];

(* ============================================================ *)
(* Rods represent apolar nematic directors: n and -n agree       *)
(* ============================================================ *)

rodSegment[pt_, n_, len_] := Line[{pt - len n, pt + len n}];

(* ============================================================ *)
(* Annulus rods panel                                            *)
(* ============================================================ *)

annulusRodPlot[alpha_?NumericQ, R_?NumericQ, dphi_?NumericQ] :=
    Module[
       {
          zp, zm, rodPts, rodList, staticLayer, rodLayer,
          rayStyle, angleArcRadius, upperAngle, lowerAngle,
          radiusLabelPos, radiusLabelOffset, angleLabelPos,
          radiusLabel, angleLabel
        },

       zp = zPlus[R, dphi];
       zm = zMinus[R, dphi];

       upperAngle = N[dphi]/2;
       lowerAngle = -N[dphi]/2;

       rayStyle = Directive[Black, AbsoluteThickness[1.0]];
       angleArcRadius = 0.58 rho;

       radiusLabelPos = 1.05 R {Cos[upperAngle], Sin[upperAngle]};
       radiusLabelOffset = 0.05 {-Sin[upperAngle], Cos[upperAngle]};
     angleLabelPos = 1. angleArcRadius {1, 0} + {-0.1, -0.18};

       radiusLabel =
    MaTeX["R", Magnification -> 1.1, "Preamble" -> maTeXPreamble];
       angleLabel =
    MaTeX["\\Delta\\phi", Magnification -> 1.1,
     "Preamble" -> maTeXPreamble];

       rodPts =
          Select[
             Flatten[

      Table[{x0, y0}, {x0, -0.9, 0.9, annulusGridSpacing}, {y0, -0.9,
        0.9, annulusGridSpacing}],
                1
              ],
             rho < Norm[#] < 1 &&
                 Norm[# - {Re[zp], Im[zp]}] > 0.08 &&
                 Norm[# - {Re[zm], Im[zm]}] > 0.08 &
           ];

       rodList =
          Table[
             Module[{n},
                n = nAnnulus[pt[[1]], pt[[2]], alpha, R, dphi];
                rodSegment[pt, n, annulusRodLength]
              ],
             {pt, rodPts}
           ];

       rodLayer = Graphics[{directorColor, directorThickness, rodList}];

       staticLayer =
          Graphics[
             {
                guideColor, guideThickness, guideDash,
                Line[{{rho, 0}, {1, 0}}],

                guideColor, guideThickness, Dashing[{}],
                Line[{{-rho, 0}, {-1, 0}}],

                innerBoundaryColor, Thick, Circle[{0, 0}, rho],
                outerBoundaryColor, Thick, Circle[{0, 0}, 1],

                rayStyle,
                Line[{{0, 0}, {Re[zp], Im[zp]}}],
                Line[{{0, 0}, {Re[zm], Im[zm]}}],

      Circle[{0, 0}, angleArcRadius, {lowerAngle, upperAngle}],

      Inset[radiusLabel, radiusLabelPos + radiusLabelOffset, Center],
                Inset[angleLabel, angleLabelPos, Center],

                defectMarker[{Re[zp], Im[zp]}, plusColor],
                defectMarker[{Re[zm], Im[zm]}, minusColor]
              }
           ];

       Show[
          rodLayer,
          staticLayer,
          PlotRange -> {{-1.05, 1.05}, {-1.05, 1.05}},
          AspectRatio -> 1,
          Frame -> False,
          Axes -> False,
          ImageSize -> annulusImageSize
        ]
     ];

(* ============================================================ *)
(* Strip rods panel                                              *)
(* ============================================================ *)

stripRodPlot[alpha_?NumericQ, R_?NumericQ, dphi_?NumericQ] :=
    Module[
       {
          nPeriods = 2,
          mMax = 3,
          vMin, vMax, u0,
          plusPts, minusPts, visiblePlusPts, visibleMinusPts,
    allDefectPts,
          xPlotRange, yPlotRange, displayScaleRatio,
          rodPts, rodList, staticLayer, rodLayer
        },

       vMin = -nPeriods Pi - N[dphi]/2 - stripEndPad;
       vMax = nPeriods Pi + N[dphi]/2 + stripEndPad;
       u0 = u0FromR[R];

       plusPts = Table[{u0, N[dphi]/2 + 2 Pi m}, {m, -mMax, mMax}];
       minusPts = Table[{u0, -N[dphi]/2 + 2 Pi m}, {m, -mMax, mMax}];

       visiblePlusPts = Select[plusPts, vMin <= #[[2]] <= vMax &];
       visibleMinusPts = Select[minusPts, vMin <= #[[2]] <= vMax &];
       allDefectPts = Join[visiblePlusPts, visibleMinusPts];

       xPlotRange = (ell + 0.15) - (-0.15);
       yPlotRange = (vMax + 0.10) - (vMin - 0.10);

       (*
   Correct rod slopes for the strip panel's compressed aspect ratio. *)
       displayScaleRatio = stripAspectRatio * xPlotRange/yPlotRange;

       rodPts =
          Flatten[
             Table[
                With[{pt = {u1, v1}},

       If[Min[Norm[pt - #] & /@ allDefectPts] > 0.12, pt, Nothing]
                 ],
                {u1, stripUSpacing/2, ell - stripUSpacing/2,
       stripUSpacing},
                {v1, vMin + stripVSpacing/2, vMax - stripVSpacing/2,
       stripVSpacing}
              ],
             1
           ];

       rodList =
          Table[
             Module[{n, nPlot},
                n = nStrip[pt[[1]], pt[[2]], alpha, R, dphi];
                nPlot = {n[[1]], n[[2]]/displayScaleRatio};
                rodSegment[pt, nPlot, stripRodLength]
              ],
             {pt, rodPts}
           ];

       rodLayer = Graphics[{directorColor, directorThickness, rodList}];

       staticLayer =
          Graphics[
             {
                guideColor, guideThickness, guideDash,
                Table[Line[{{0, 2 Pi m}, {ell, 2 Pi m}}], {m, -2, 2}],

                guideColor, guideThickness, Dashing[{}],

      Table[Line[{{0, (2 m + 1) Pi}, {ell, (2 m + 1) Pi}}], {m, -2, 1}],

                innerBoundaryColor, Thick, Line[{{0, vMin}, {0, vMax}}],
                outerBoundaryColor, Thick,
      Line[{{ell, vMin}, {ell, vMax}}],


      Table[defectMarker[pt, plusColor], {pt, visiblePlusPts}],

      Table[defectMarker[pt, minusColor], {pt, visibleMinusPts}]
              }
           ];

       Show[
          rodLayer,
          staticLayer,

    PlotRange -> {{-0.15, ell + 0.15}, {vMin - 0.10, vMax + 0.10}},
          AspectRatio -> stripAspectRatio,
          Frame -> False,
          Axes -> False,
          ImageSize -> stripImageWidth
        ]
     ];

(* ============================================================ *)
(* Upper-half-plane rods panel                                   *)
(* ============================================================ *)

uhpRodPlot[alpha_?NumericQ, R_?NumericQ, dphi_?NumericQ] :=
    Module[
       {zetap, zetam, zetapImg, zetamImg, yBelow, rodPts, rodList,
    staticLayer, rodLayer},

       zetap = zetaPlus[R, dphi];
       zetam = zetaMinus[R, dphi];
       zetapImg = zetaPlusImage[R, dphi];
       zetamImg = zetaMinusImage[R, dphi];

       yBelow = uhpBelowAxisFactor uhpYMax;

       rodPts =
          Flatten[
             Table[
                With[{pt = {x0, y0}},
                   If[
                      Norm[pt - {Re[zetap], Im[zetap]}] > 0.20 &&
                       Norm[pt - {Re[zetam], Im[zetam]}] > 0.20,
                      pt,
                      Nothing
                    ]
                 ],
                {x0, -uhpXMax + uhpXSpacing/2,
       uhpXMax - uhpXSpacing/2, uhpXSpacing},
                {y0, uhpYSpacing/2, uhpYMax - uhpYSpacing/2,
       uhpYSpacing}
              ],
             1
           ];

       rodList =
          Table[
             Module[{n},
                n = nUHP[pt[[1]], pt[[2]], alpha, R, dphi];
                rodSegment[pt, n, uhpRodLength]
              ],
             {pt, rodPts}
           ];

       rodLayer = Graphics[{directorColor, directorThickness, rodList}];

       staticLayer =
          Graphics[
             {
                guideColor, guideThickness, Dashed,
                Circle[{0, 0}, 1, {0, Pi}],

                Dashing[{}], innerBoundaryColor, Thick,
      Line[{{-uhpXMax, 0}, {0, 0}}],
                outerBoundaryColor, Thick, Line[{{0, 0}, {uhpXMax, 0}}],

                clusterMarker[{0, 0}],
                defectMarker[{Re[zetap], Im[zetap]}, plusColor],
                defectMarker[{Re[zetam], Im[zetam]}, minusColor],

      defectMarker[{Re[zetapImg], Im[zetapImg]}, plusColor,
       imageOpacity],

      defectMarker[{Re[zetamImg], Im[zetamImg]}, minusColor,
       imageOpacity]
              }
           ];

       Show[
          rodLayer,
          staticLayer,
          PlotRange -> {{-uhpXMax, uhpXMax}, {-yBelow, uhpYMax}},
          AspectRatio -> (uhpYMax + yBelow)/(2 uhpXMax),
          Frame -> False,
          Axes -> False,
          ImageSize -> uhpImageWidth
        ]
     ];

(* ============================================================ *)
(* Build the three panels                                       *)
(* ============================================================ *)

annulusPanel = annulusRodPlot[alpha, R, dphi];
stripPanel = stripRodPlot[alpha, R, dphi];
uhpPanel = uhpRodPlot[alpha, R, dphi];

baseFigure =
    Row[
       {
          annulusPanel,
          Spacer[gapAnnulusStrip],
          stripPanel,
          Spacer[gapStripUHP],
          uhpPanel
        },
       Alignment -> Center
     ];

figureWidth =
  annulusImageSize + gapAnnulusStrip + stripImageWidth + gapStripUHP +
   uhpImageWidth;

figureHeight =
    Max[
       annulusImageSize,
       stripImageWidth stripAspectRatio,
       uhpImageWidth (uhpYMax + uhpBelowAxisFactor uhpYMax)/(2 uhpXMax)
     ];

(* ============================================================ *)
(* Annotation layout: reference pixel coordinates scaled to the  *)
(* actual composite height; fixed horizontal positions retained. *)
(* ============================================================ *)

referenceStripAspectRatio = 3.5;
referenceFigureHeight = stripImageWidth referenceStripAspectRatio;
annotationYScale = figureHeight/referenceFigureHeight;

scaleY[y_] := annotationYScale y;
scalePoint[{x_, y_}] := {x, scaleY[y]};

zLabel = MaTeX["z", Magnification -> planeMagnification,
   "Preamble" -> maTeXPreamble];
ZLabel = MaTeX["Z", Magnification -> planeMagnification,
   "Preamble" -> maTeXPreamble];
zetaLabel =
  MaTeX["\\zeta", Magnification -> planeMagnification,
   "Preamble" -> maTeXPreamble];

mapEquation1 =
  MaTeX["Z=\\log(z/\\rho)", Magnification -> mapMagnification,
   "Preamble" -> maTeXPreamble];
mapEquation2 =
  MaTeX["\\zeta=-\\exp\\!\\left[-\\frac{\\mathrm{i}\\pi \
Z}{\\log(1/\\rho)}\\right]",
     Magnification -> mapMagnification,
     "Preamble" -> maTeXPreamble
   ];

zPlanePosReference = {190, 360};
ZPlanePosReference = {380, 490};
zetaPlanePosReference = {575, 420};
mapEquation1PosReference = {170, 410};
mapEquation2PosReference = {480, 450};

arrow1PtsReference = {{205, 380}, {217, 405}, {233, 415}, {247, 424}};
arrow2PtsReference = {{385, 480}, {385, 462}, {385, 440}, {420, 415}};

zPlanePos = scalePoint[zPlanePosReference];
ZPlanePos = scalePoint[ZPlanePosReference];
zetaPlanePos = scalePoint[zetaPlanePosReference];
mapEquation1Pos = scalePoint[mapEquation1PosReference];
mapEquation2Pos = scalePoint[mapEquation2PosReference];
arrow1Pts = scalePoint /@ arrow1PtsReference;
arrow2Pts = scalePoint /@ arrow2PtsReference;

planeSX = 17;
planeSY = 17;

mapOverlay =
    Graphics[
       {

    planeMarker[zPlanePos, zLabel, planeSX, planeSY, {0.52, 0.61}],

    planeMarker[ZPlanePos, ZLabel, planeSX, planeSY, {0.52, 0.61}],

    planeMarker[zetaPlanePos, zetaLabel, planeSX,
     planeSY, {0.54, 0.70}],

          Inset[mapEquation1, mapEquation1Pos, Center],
          Inset[mapEquation2, mapEquation2Pos, Center],

          Black,
          AbsoluteThickness[1.25],
          Arrowheads[0.025],
          Arrow[BezierCurve[arrow1Pts]],
          Arrow[BezierCurve[arrow2Pts]]
        },
       PlotRange -> {{0, figureWidth}, {0, figureHeight}},
       PlotRangePadding -> None,
       ImagePadding -> 0,
       Background -> None,
       ImageSize -> {figureWidth, figureHeight}
     ];

conformalMapFigure =
  Overlay[{baseFigure, mapOverlay}, Alignment -> Center];

(* ============================================================ *)
(* Export                                                       *)
(* ============================================================ *)

Export["Fig2.pdf", conformalMapFigure, "PDF",
  ImageResolution -> 150]

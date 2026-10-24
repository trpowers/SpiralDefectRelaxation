(* Figures 7, 8, and 9: Circular-annulus defect energetics.
   A single script retains the common calculations from the source notebook.

   Fig7.pdf: optimal radial position (formerly R-star-v-Dphi.pdf).
   Fig8.pdf: relative free energy and barriers (FtotWithBarriers.pdf).
   Fig9.pdf: critical angular separation (BarrierPosition.pdf).

   Requires Mathematica and MaTeX with its LaTeX dependencies.
   Exports all three PDFs to the current working directory, Directory[].
   As in the original notebook, the following Clear removes Global` values.
*)

Needs["MaTeX`"];
Clear["Global`*"];

(* 1. Shared parameters and plot styles.
   L is ell = Log[1/rho] in the manuscript. The energy and barrier plots
   use L = 1; Fig7 varies L over Lvals. Energies below are divided by K. *)

cbPalette = {RGBColor[0, 0.447, 0.698],(*blue*)
   RGBColor[0.835, 0.369, 0],(*orange*)
   RGBColor[0, 0.620, 0.451],(*bluish green*)
   RGBColor[0.800, 0.475, 0.655],(*reddish purple*)
   RGBColor[0.55, 0.35, 0.00],(*brown*)
   RGBColor[0.337, 0.706, 0.914]   (*sky blue*)};

cbStyles = Directive[#, Thick] & /@ cbPalette;

L = 1;
epsOverRho = 0.02;
cCore = Pi/8;
alphaVals = Range[Pi/6, Pi, Pi/6];
nSGrid = 80;
dphiMin = 3 epsOverRho;

(* Geometry: s = 1 + Log[R]/L runs from 0 to 1 across the annulus.
   epsPhys is the physical core cutoff, with epsPhys/rho = epsOverRho.
   The radial search stays two cutoff radii away from each boundary. *)

rho[L_] := Exp[-L];
q[L_] := Exp[-Pi^2/L];

epsPhys[L_] := epsOverRho rho[L];

RofS[s_, L_] := Exp[L (s - 1)];
sOfR[R_, L_] := 1 + Log[R]/L;

omegaOfS[s_, L_] := Exp[2 Pi I Log[RofS[s, L]]/L];

Pcrowdy[lambda_?NumericQ, qq_?NumericQ] := 
  QPochhammer[lambda, qq^2] QPochhammer[qq^2/lambda, qq^2];

sMin[L_] := sOfR[rho[L] + 2 epsPhys[L], L];
sMax[L_] := sOfR[1 - 2 epsPhys[L], L];

(* Exact Crowdy product, pair energy, and radial minimization.
   Xfun is the manuscript's X. Keep the original finite-core expression,
   overlap guard, grid search, local refinement, and machine precision. *)

Xfun[x_?NumericQ, om_?NumericQ, qq_?NumericQ] := 
  Pcrowdy[x, qq]*Pcrowdy[om*x, qq]*Pcrowdy[1/x, qq]*Pcrowdy[om/x, qq];

fPairExactOverK[s_?NumericQ, dphi_?NumericQ, L_?NumericQ] := 
 Module[{R, qq, eps, delta, a, b, omega}, R = RofS[s, L];
  qq = q[L];
  eps = epsPhys[L];
  delta = eps/R;
  omega = omegaOfS[s, L];
  If[dphi <= 2 delta, Return[10^20]];
  a = Exp[Pi delta/L];
  b = Exp[-Pi (dphi - delta)/L];
  -(Pi/4) Log[Abs[Xfun[a, omega, qq]/Xfun[b, omega, qq]]]];

sStarExact[dphi_?NumericQ, L_?NumericQ] := 
  Module[{sGrid, vals, goodVals, pos, s0, sLeft, sRight, s, res}, 
   sGrid = Subdivide[sMin[L], sMax[L], nSGrid];
   vals = Table[{s, fPairExactOverK[s, dphi, L]}, {s, sGrid}];
   goodVals = Select[vals, NumericQ[#[[2]]] && #[[2]] < 10^10 &];
   pos = First@Ordering[goodVals[[All, 2]], 1];
   s0 = goodVals[[pos, 1]];
   sLeft = If[pos == 1, goodVals[[pos, 1]], goodVals[[pos - 1, 1]]];
   sRight = 
    If[pos == Length[goodVals], goodVals[[pos, 1]], 
     goodVals[[pos + 1, 1]]];
   res = 
    Quiet@FindMinimum[
      fPairExactOverK[s, dphi, L], {s, s0, sLeft, sRight}, 
      WorkingPrecision -> MachinePrecision];
   {s /. Last[res], First[res]}];

(* 2. Shared L = 1 data for Figures 8 and 9.
   Preserve the original 281-point separation grid, clustered near 2 Pi.
   Columns: {DeltaPhi/Pi, DeltaPhi, sStar, RStar, Fpair/K}. *)

dphiValsEnergy = 
  DeleteDuplicates[
   Sort[Join[Subdivide[dphiMin, 1.8 Pi, 120], 
     2 Pi - (Exp /@ 
        Subdivide[Log[2 Pi - 1.8 Pi], Log[3 epsOverRho], 160])]]];

rstarData = Table[
  Module[{ss, ff},
    {ss, ff} = sStarExact[dphi, L];
    {dphi/Pi, dphi, ss, RofS[ss, L], ff}
  ],
  {dphi, dphiValsEnergy}
];

(* Relative energy: DeltaF/K = 2 cCore + Fpair/K - Pi alpha DeltaPhi/L. *)

deltaFtotOverK[dphi_, fpair_, alpha_] := 
  2 cCore + fpair - (Pi alpha/L) dphi;

deltaFCurve[alpha_] := 
  rstarData /. {x_, dphi_, ss_, rr_, fpair_} :> {x, 
     deltaFtotOverK[dphi, fpair, alpha]};

(* 3. Figure 7: optimal radial positions for four annulus widths.
   The source notebook used a separate 241-point grid here. Keeping it
   separate avoids silently changing either this plot or Figures 8 and 9. *)

Lvals = {1/2, 1, 2, 3};

dphiValsRadius = 
  DeleteDuplicates@
   Sort@Join[Subdivide[dphiMin, 1.8 Pi, 120], 
     2 Pi - (Exp /@ 
        Subdivide[Log[2 Pi - 1.8 Pi], Log[3 epsOverRho], 120])];

(* For each L, columns are {DeltaPhi/Pi, sStar, RStar, Fpair/K}. *)

dataForL[L_?NumericQ] := 
  Table[Module[{ss, ff}, {ss, ff} = 
     sStarExact[dphi, L];{dphi/Pi, ss, 
     RofS[ss, L], ff}], {dphi, dphiValsRadius}];

rawData = Table[dataForL[L], {L, Lvals}];

sStarData = rawData[[All, All, {1, 2}]];

plots = ListLinePlot[sStarData, 
   PlotStyle -> cbStyles[[;; Length[Lvals]]], InterpolationOrder -> 1];

pltRstarExact = 
  Show[plots, Frame -> True, 
   FrameLabel -> {MaTeX["\\Delta\\phi/\\pi", Magnification -> 1.6], 
     MaTeX["1-\\log R_{*}/\\log \\rho", Magnification -> 1.6]}, 
   LabelStyle -> {FontFamily -> "Times New Roman", FontSize -> 16}, 
   PlotRange -> {{0, 2}, {0, 1/2}}, PlotRangePadding -> None, 
   AspectRatio -> 0.8, ImageSize -> 500];

(* 4. Barrier calculation shared by Figures 8 and 9.
   Retain the original interpolation and search interval. The search is
   for a maximum on the first half of the path, not an energy zero. *)

criticalSeparation[alpha_?NumericQ] := 
  Module[{curve, interp, xMin, xMax, res, xc, fc}, 
   curve = deltaFCurve[alpha];
   interp = Interpolation[curve];
   xMin = Max[First[curve[[All, 1]]], 0.03];
   xMax = 0.98;
   res = Quiet@NMaximize[{interp[x], xMin <= x <= xMax}, x];
   xc = x /. Last[res];
   fc = First[res];
   {xc, fc}];

(* Figure 9 starts just above alphaCrit, as in the source notebook.
   Columns: {alpha/Pi, DeltaPhiCritical/Pi, DeltaFCritical/K}. *)

alphaCrit = Pi/2;

alphaValsCritical = 
  DeleteDuplicates@
   N@Join[Range[alphaCrit + Pi/1000, 0.60 Pi, Pi/1000], 
     Range[0.60 Pi, Pi, Pi/200]];

criticalData = 
  N@Table[Module[{xc, fc}, {xc, fc} = criticalSeparation[alpha];
     {alpha/Pi, xc, fc}], {alpha, alphaValsCritical}];

pltCriticalSeparation = 
  Show[ListLinePlot[criticalData[[All, {1, 2}]], Frame -> True, 
    FrameLabel -> {MaTeX["\\alpha/\\pi", Magnification -> 1.8], 
      MaTeX["\\Delta\\phi_{\\mathrm{c}}/\\pi", Magnification -> 1.8]},
     PlotStyle -> Directive[cbPalette[[1]], Thick], 
    LabelStyle -> {FontFamily -> "Times New Roman", FontSize -> 18}, 
    PlotRange -> {{alphaCrit/Pi, 1}, {0, 1/2}}, 
    PlotRangePadding -> Scaled[0.03], AspectRatio -> 0.8, 
    ImageSize -> 500], 
   Graphics[{Black, Dashed, Thin, 
     Line[{{alphaCrit/Pi, 0}, {alphaCrit/Pi, 1/2}}]}]];

(* 5. Figure 8: energy curves and their barrier markers.
   The original fit over 1.05 <= DeltaPhi/Pi <= 1.60 is needed for the
   marker-selection rule. Its diagnostic plot and printed output are not. *)

gCurve = 
  rstarData /. {x_, dphi_, ss_, rr_, fpair_} :> {x, 2 cCore + fpair};

(*Fit the nearly straight part of G*)
xTailMin = 1.05;
xTailMax = 1.60;

tailData = Select[gCurve, xTailMin <= #[[1]] <= xTailMax &];

tailFit = LinearModelFit[tailData, x, x];

gSlopeX = tailFit["BestFitParameters"][[2]];

(* alphaNuc is the original tail-slope estimate. It is retained only for
   selecting the curves that receive barrier dots, not as a replacement
   for the manuscript's symmetry threshold alphaCrit = Pi/2. *)
alphaNuc = L gSlopeX/Pi^2;

alphaValsDisplay = alphaVals;

dataFDisplay = deltaFCurve /@ alphaValsDisplay;

pltFtotExact = 
  Show[ListLinePlot[dataFDisplay, 
    PlotStyle -> cbStyles[[;; Length[alphaValsDisplay]]], 
    InterpolationOrder -> 1], Frame -> True, 
   FrameLabel -> {MaTeX["\\Delta\\phi/\\pi", Magnification -> 1.8], 
     MaTeX["\\Delta F/K", Magnification -> 1.8]}, 
   LabelStyle -> {FontFamily -> "Times New Roman", FontSize -> 18}, 
   PlotRange -> {{0, 2}, All}, PlotRangePadding -> Scaled[0.03], 
   AspectRatio -> 0.8, ImageSize -> 500];

alphaValsForCriticalDots = Select[alphaValsDisplay, # > alphaNuc &];

criticalPointsForAlphaVals = 
  Table[Module[{xc, fc}, {xc, fc} = criticalSeparation[alpha];
    {xc, fc}], {alpha, alphaValsForCriticalDots}];

pltFtotWithCriticalDots = 
  Show[pltFtotExact, 
   ListPlot[criticalPointsForAlphaVals, 
    PlotStyle -> {Black, PointSize[0.010]}]];

(* 6. Export only the three figures used in the manuscript. *)
Export["Fig7.pdf", pltRstarExact, "PDF", ImageResolution -> 150];
Export["Fig8.pdf", pltFtotWithCriticalDots, "PDF", ImageResolution -> 150];
Export["Fig9.pdf", pltCriticalSeparation, "PDF", ImageResolution -> 150];

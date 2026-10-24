(* Figure 12: Free energy and barrier markers for a needle-like inclusion.
   Original figure: FtotWithBarriersNeedle.pdf.

   Requires Mathematica and MaTeX with its LaTeX dependencies.
   Run this file as a package, or load it with Get from a notebook.
   Exports Fig12.pdf to the current working directory, Directory[].
   As in the source notebook, Clear below removes existing Global` values.

   Retains the original energy formula, separation grid, minimization
   settings, interpolation, and styling. The very-small-separation samples
   approach core contact and should not be interpreted as resolved
   nucleation energies. Only the manuscript's energy plot is exported.
*)

Needs["MaTeX`"];
Clear["Global`*"];

(* 1. Plot styles. Colors follow the six alpha values in increasing order. *)

cbPalette = {RGBColor[0, 0.447, 0.698], RGBColor[0.835, 0.369, 0], 
   RGBColor[0, 0.620, 0.451], RGBColor[0.800, 0.475, 0.655], 
   RGBColor[0.55, 0.35, 0.00], RGBColor[0.337, 0.706, 0.914]};

cbStyles = Directive[#, Thick] & /@ cbPalette;

(* 2. Confocal elliptical annulus: z = c Exp[I alpha] Cosh[u + I v].
   c is the center-to-focus distance; L is ell = uOuter - uInner.
   epsPhys is a fixed physical core radius, not an auxiliary-annulus radius.
   Each defect has Fcore/K = cCore. All energies below are divided by K. *)

c = 1;
uInner = 0.1;
uOuter = 2.0;

L = uOuter - uInner;
q = Exp[-Pi^2/L];

epsOverInnerMinor = 0.02;
cCore = Pi/8;

bInner = c Sinh[uInner];
epsPhys = epsOverInnerMinor bInner;

nSGrid = 80;

(* s = (u - uInner)/L. R = Exp[u - uOuter] is the radius in the auxiliary
   circular annulus, not a physical distance in the elliptical geometry.
   hScale = |dz/dW| converts coordinate increments into physical lengths.
   The modulus-one factor omega equals Exp[2 Pi I (u - uInner)/L]. *)

hScale[u_, v_] := c/Sqrt[2] Sqrt[Cosh[2 u] - Cos[2 v]];

RofS[s_] := Exp[L (s - 1)];
uOfS[s_] := uInner + L s;

omegaOfS[s_] := Exp[2 Pi I Log[RofS[s]]/L];

(* 3. Exact prime-function product, with the manuscript's nome q.
   P(lambda,q) = (lambda;q^2)_infinity (q^2/lambda;q^2)_infinity.
   Xfun is the same four-factor product as in the circular calculation.
   "Exact" refers to P; the finite-core treatment remains the manuscript's
   small-core approximation. *)

Pcrowdy[lambda_?NumericQ, qq_?NumericQ] := 
  QPochhammer[lambda, qq^2] QPochhammer[qq^2/lambda, qq^2];

Xfun[x_?NumericQ, om_?NumericQ, qq_?NumericQ] := 
  Pcrowdy[x, qq]*Pcrowdy[om*x, qq]*Pcrowdy[1/x, qq]*Pcrowdy[om/x, qq];

(* 4. Original 281-point separation grid, clustered near 2 Pi.
   dv = vPlus - vMinus, with vPlus = dv/2 and vMinus = -dv/2.
   It is an elliptic-coordinate separation, not the lab polar-angle gap.
   The pair starts near one tip (dv -> 0) and meets near the other (dv -> 2 Pi).
   Both coalescence endpoints are excluded. The source grid is preserved;
   core clearance must be rechecked if the geometry or sampling is changed. *)

dvMin = 3 epsPhys/hScale[uInner + 0.1 L, 0];

dvVals = 
  DeleteDuplicates[
   Sort[Join[Subdivide[dvMin, 1.8 Pi, 120], 
     2 Pi - (Exp /@ Subdivide[Log[2 Pi - 1.8 Pi], Log[3 dvMin], 160])]]];

(* 5. Pair energy for a trial s.
   The common local angular core cutoff is delta = epsPhys/hScale[u,dv/2].
   h has the same value at the two symmetric defect positions.
   a = Exp[Pi delta/L], b = Exp[-Pi (dv - delta)/L].
   This implements Fpair/K = -(Pi/4) Log[Abs[X(a)/X(b)]]. *)

fPairNeedleOverK[s_?NumericQ, dv_?NumericQ] := 
 Module[{uStar, vPlus, h, delta, a, b, omega, qq},
  uStar = uOfS[s];
  vPlus = dv/2;
  h = hScale[uStar, vPlus];
  delta = epsPhys/h;
  qq = q;
  omega = omegaOfS[s];
  (*Original local core-overlap and boundary-clearance guards.*)
  If[dv <= 2 delta, Return[10^20]];
  If[uStar - uInner <= 2 delta, Return[10^20]];
  If[uOuter - uStar <= 2 delta, Return[10^20]];
  a = Exp[Pi delta/L];
  b = Exp[-Pi (dv - delta)/L];
  -(Pi/4) Log[Abs[Xfun[a, omega, qq]/Xfun[b, omega, qq]]]];

(* 6. Corrected coarse-grid search followed by local refinement.
   It uses the ORIGINAL neighboring grid points even when one was rejected.
   At a rejected neighbor, bisection locates the existing admissibility
   boundary from within the allowed region (tolerance edgeTol = 10^-10 in s).
   The energies at the admissible search endpoints are also compared with
   the local FindMinimum result. No physical cutoff is changed.
   This is the corrected minimizer verified with the working notebook:
   first saved sample s* = 0.0497676 and Fpair/K approximately 2.3*10^-9.
   Near this sample the defect cores almost touch: the point-defect model
   does not quantitatively resolve nucleation or annihilation. *)

sMin = 0.02;
sMax = 0.98;

sStarNeedle[dv_?NumericQ] := 
 Module[{sGrid, vals, goodVals, pos, s0, sLeft, sRight, s, res, 
   gridPos, allowedEdge, sStart, sFit, candidates, edgeTol = 10^-10}, 
  sGrid = Subdivide[sMin, sMax, nSGrid];
  vals = Table[{s, fPairNeedleOverK[s, dv]}, {s, sGrid}];
  goodVals = Select[vals, NumericQ[#[[2]]] && #[[2]] < 10^10 &];
  If[Length[goodVals] < 3, Return[{Missing[], 10^20}]];
  pos = First@Ordering[goodVals[[All, 2]], 1];
  s0 = goodVals[[pos, 1]];
  
  gridPos = First@FirstPosition[sGrid, s0];
  sLeft = sGrid[[Max[1, gridPos - 1]]];
  sRight = sGrid[[Min[Length[sGrid], gridPos + 1]]];
  
  allowedEdge[trial_, inside_] := 
   Module[{bad = trial, good = inside, mid}, 
    If[TrueQ[fPairNeedleOverK[bad, dv] < 10^10], Return[bad]];
    While[Abs[good - bad] > edgeTol, mid = (bad + good)/2;
     If[TrueQ[fPairNeedleOverK[mid, dv] < 10^10], good = mid, 
      bad = mid]];
    good];
  sLeft = allowedEdge[sLeft, s0];
  sRight = allowedEdge[sRight, s0];
  
  candidates = {goodVals[[pos]], {sLeft, 
     fPairNeedleOverK[sLeft, dv]}, {sRight, 
     fPairNeedleOverK[sRight, dv]}};
  sStart = If[sLeft < s0 < sRight, s0, (sLeft + sRight)/2];
  res = Quiet@
    FindMinimum[fPairNeedleOverK[s, dv], {s, sStart, sLeft, sRight}, 
     WorkingPrecision -> MachinePrecision];
  If[MatchQ[res, {_?NumericQ, {__Rule}}], sFit = s /. Last[res];
   If[TrueQ[sLeft <= sFit <= sRight], 
    AppendTo[candidates, {sFit, fPairNeedleOverK[sFit, dv]}]]];
  candidates = 
   Select[candidates, NumericQ[#[[2]]] && TrueQ[#[[2]] < 10^10] &];
  candidates[[First@Ordering[candidates[[All, 2]], 1]]]];

(* 7. Compute the shared optimized-energy data once.
   Columns: {dv/Pi, dv, sStar, RStar, uStar, Fpair/K}.
   RStar and uStar are retained to make the intermediate data inspectable. *)

rstarData = Table[
  Module[{ss, ff},
    {ss, ff} = sStarNeedle[dv];
    If[ss === Missing[], Nothing,
      {dv/Pi, dv, ss, RofS[ss], uOfS[ss], ff}]
  ],
  {dv, dvVals}
];

(* Stop rather than export an incomplete or invalid dataset.
   No diagnostic messages are printed on a successful run. *)
fig12::data = "The radial calculation did not return a complete numerical dataset; no PDF was exported.";
If[Length[rstarData] != Length[dvVals] ||
    !MatrixQ[rstarData, NumericQ],
  Message[fig12::data];
  Abort[]
];
If[Max[rstarData[[All, 6]]] >= 10^10,
  Message[fig12::data];
  Abort[]
];

(* 8. Relative energy: subtract the defect-free spiral energy.
   DeltaF/K = 2 cCore + Fpair/K - Pi alpha dv/L.
   The drawn curves retain the original linear interpolation and styling. *)

deltaFtotNeedleOverK[dv_, fpair_, alpha_] := 
  2 cCore + fpair - (Pi alpha/L) dv;

deltaFCurve[alpha_] := 
  rstarData /. {x_, dv_, ss_, rr_, uu_, fpair_} :> {x, 
     deltaFtotNeedleOverK[dv, fpair, alpha]};

alphaVals = Range[Pi/6, Pi, Pi/6];

dataF = deltaFCurve /@ alphaVals;

pltFtotNeedle = 
  Show[ListLinePlot[dataF, PlotStyle -> cbStyles[[;; Length[alphaVals]]],
     InterpolationOrder -> 1], Frame -> True, 
   FrameLabel -> {MaTeX["\\Delta v/\\pi", Magnification -> 1.8], 
     MaTeX["\\Delta F/K", Magnification -> 1.8]}, 
   LabelStyle -> {FontFamily -> "Times New Roman", FontSize -> 18}, 
   PlotRange -> {{0, 2}, All}, PlotRangePadding -> Scaled[0.03], 
   AspectRatio -> 0.8, ImageSize -> 500];

(* 9. Barrier markers only for the displayed alpha values.
   alphaCrit is supplied by the manuscript's symmetry argument, not fitted.
   For alpha > Pi/2, retain the original cubic Interpolation and maximization
   interval in x = dv/Pi. This locates the first barrier for these curves.
   At alpha = Pi/2, the original code places the marker at x = 1 by symmetry;
   it does not obtain that point from the restricted maximization below.
   The finite-cutoff formula may have small higher-order symmetry errors. *)

alphaCrit = Pi/2;

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

alphaValsForCriticalDots = Select[alphaVals, # >= alphaCrit &];

criticalPointsForAlphaVals = 
  Table[If[Chop[alpha - alphaCrit] == 0, 
    With[{interp = Interpolation[deltaFCurve[alpha]]}, {1, interp[1]}],
     criticalSeparation[alpha]], {alpha, alphaValsForCriticalDots}];

pltFtotNeedleWithCriticalDots = 
  Show[pltFtotNeedle, 
   ListPlot[criticalPointsForAlphaVals, 
    PlotStyle -> {Black, PointSize[0.010]}]];

(* 10. Export only the manuscript figure, to Directory[]. *)

Export["Fig12.pdf", pltFtotNeedleWithCriticalDots, "PDF",
  ImageResolution -> 150];

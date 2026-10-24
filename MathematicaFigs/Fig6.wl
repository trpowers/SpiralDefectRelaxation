(* Figure 6: Integration contours around a defect pair in an annulus.
   Reproduces the original graphic using the original approximate
   prime-function representation and the same rod and contour geometry.

   Dependencies: Mathematica, MaTeX` (and its LaTeX installation).
   Output: Fig6.pdf in Directory[].
*)

Needs["MaTeX`"]

ClearAll[Omega, twoIOmega, rawTheta, nfieldDef, makeDefectPlot, Zmap, \
Papprox]

(* Geometry and director-field parameters. *)
rho = 1/E;
a = rho;
b = 1;
H = b;

L = Log[1/rho];
q = Exp[-Pi^2/L];

rPlus = 0.7;
rMinus = 0.7;

diskRad = 0.05;
epsCore = 0.065;
cutSep = 0.018;

(* Offset of the magenta integration contours from the boundaries. *)
outerContourSep = 0.018;

(* Equivalent to the manuscript map Z = Log[z/rho],
   zeta = -Exp[-I Pi Z/L], because Log[rho] = -L. *)
Zmap[z_] := Exp[-I Pi Log[z]/L];

(* Leading Crowdy product and the first correction on each side.
   For L = 1, q^2 = Exp[-2 Pi^2] is approximately 2.7 x 10^-9.
   This approximation is part of the original director-rod illustration;
   it should not be substituted for the exact function in energy work. *)
Papprox[lambda_] := (1 - lambda) (1 - q^2 lambda) (1 - q^2/lambda);

makeDefectPlot[phi_, alpha_ : 0] := 
 Module[{phiPlus, phiMinus, zPlus, zMinus, zetaPlus, zetaMinus, 
   twoIOmega, Omega, rawTheta, thetaPlot, nfieldDef, rCutPlus, 
   rCutMinus, trimPlus, trimMinus, gapPlus, gapMinus, gammaPlusAngle, 
   rodPts, rodLength, rodList, rOuterContour, rInnerContour, 
   boundaryContours}, phiPlus = phi;
  phiMinus = -phi;
  zPlus = rPlus Exp[I phiPlus];
  zMinus = rMinus Exp[I phiMinus];
  zetaPlus = Zmap[zPlus];
  zetaMinus = Zmap[zMinus];
  twoIOmega[z_] := 
   Module[{zeta, lambdaPlus, lambdaMinus, lambdaTildePlus, 
     lambdaTildeMinus}, zeta = Zmap[z];
    lambdaPlus = zetaPlus/zeta;
    lambdaMinus = zetaMinus/zeta;
    lambdaTildePlus = Conjugate[zetaPlus]/zeta;
    lambdaTildeMinus = Conjugate[zetaMinus]/zeta;
    Log[(Papprox[lambdaPlus] Papprox[lambdaTildePlus])/(Papprox[
         lambdaMinus] Papprox[lambdaTildeMinus])]];
  Omega[z_] := -I/2 twoIOmega[z];
  rawTheta[x_, y_] := Re[Omega[x + I y]];
  thetaPlot[x_, y_] := 
   rawTheta[x, y] + alpha Log[Sqrt[x^2 + y^2]]/Log[rho];
  nfieldDef[x_, y_] := {Cos[thetaPlot[x, y]], Sin[thetaPlot[x, y]]};
  rCutPlus = rPlus + cutSep;
  rCutMinus = rPlus - cutSep;
  trimPlus = 
   ArcCos[(rPlus^2 + rCutPlus^2 - epsCore^2)/(2 rPlus rCutPlus)];
  trimMinus = 
   ArcCos[(rPlus^2 + rCutMinus^2 - epsCore^2)/(2 rPlus rCutMinus)];
  gapPlus = ArcSin[cutSep/epsCore];
  gapMinus = ArcSin[cutSep/epsCore];
  gammaPlusAngle = 0.12;
  (* Full inner and outer boundary contours, inside the annulus. *)
  rOuterContour = b - outerContourSep;
  rInnerContour = a + outerContourSep;
  boundaryContours = {Circle[{0, 0}, rOuterContour], 
    Circle[{0, 0}, rInnerContour]};
  rodLength = 0.045;
  rodPts = 
   Select[Flatten[
     Table[{x0, y0}, {x0, -0.9, 0.9, 0.18}, {y0, -0.9, 0.9, 0.18}], 1],
     a < Norm[#] < b &];
  rodList = 
   Table[Module[{x0 = pt[[1]], y0 = pt[[2]], n}, n = nfieldDef[x0, y0];
     Line[{pt - rodLength n, pt + rodLength n}]], {pt, rodPts}];
  Graphics[{RGBColor[0.29, 0.46, 0.6], AbsoluteThickness[1.0], 
    rodList, Black, Thick, Circle[{0, 0}, a], Circle[{0, 0}, b], Blue,
     Disk[{Re[zPlus], Im[zPlus]}, diskRad], Red, 
    Disk[{Re[zMinus], Im[zMinus]}, diskRad], Magenta, Thick, 
    boundaryContours, 
    Circle[{Re[zPlus], Im[zPlus]}, 
     epsCore, {phiPlus - Pi/2 + gapPlus, phiPlus + 3 Pi/2 - gapPlus}],
     Circle[{Re[zMinus], Im[zMinus]}, 
     epsCore, {phiMinus + Pi/2 + gapMinus, 
      phiMinus + 5 Pi/2 - gapMinus}], 
    Circle[{0, 0}, rCutPlus, {-phi + trimPlus, phi - trimPlus}], 
    Circle[{0, 0}, rCutMinus, {-phi + trimMinus, phi - trimMinus}], 
    Text[MaTeX["\\Gamma_{-}", 
      Magnification -> 
       1.25], (rCutPlus + 0.15) {Cos[0.9 gammaPlusAngle], 
       Sin[0.9 gammaPlusAngle]}], 
    Text[MaTeX["\\Gamma_{+}", 
      Magnification -> 
       1.25], (rCutMinus - 0.1) {Cos[1.1 gammaPlusAngle], 
       Sin[1.1 gammaPlusAngle]}]}, Frame -> False, Axes -> False, 
   AspectRatio -> 1, PlotRange -> {{-H, H}, {-H, H}}, 
   PlotRangePadding -> Scaled[0.015], ImagePadding -> 3, 
   ImageSize -> 220]]

(* The manuscript illustration has phi = Pi/3 and alpha = 3 Pi/4. *)
pltC = makeDefectPlot[Pi/3, 3 Pi/4];

Export["Fig6.pdf", pltC, "PDF", ImageResolution -> 150];

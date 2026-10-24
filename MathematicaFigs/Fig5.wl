(* Figure 5: Director textures for a pair of +/-1/2 defects in an annulus.

   Geometry: rho < |z| < 1, with rho = Exp[-1].
   Defects: z+ = R Exp[I phi], z- = R Exp[-I phi]; DeltaPhi = 2 phi.
   Rows: alpha = Pi/4, Pi/2, 3 Pi/4.
   Columns: phi = Pi/12, Pi/3, 2 Pi/3, 11 Pi/12.

   Manuscript map: Z = Log[z/rho], zeta = -Exp[-I Pi Z/L],
   with L = Log[1/rho]. Because Log[rho] = -L,
       -Exp[-I Pi Log[z/rho]/L] = Exp[-I Pi Log[z]/L].
   We use the equivalent expression on the right below.

   The director angle is the sum of the defect field (the principal
   logarithm of the prime-function ratio) and the spiral field.

   Output: Fig5.pdf in the current working directory.
*)

ClearAll[rho, a, b, H, L, q, rPlus, rMinus, Zmap, Pcrowdy,
  makeDefectPlot, phiVals, alphaVals, pltGridAlpha, x, y, z, phi, alpha];

(* Geometry and prime-function parameter. *)
rho = 1/E;
a = rho;
b = 1;
H = b;
L = Log[1/rho];
q = Exp[-Pi^2/L];

rPlus = 0.7;
rMinus = 0.7;

Zmap[z_] := Exp[-I Pi Log[z]/L];
Pcrowdy[lambda_] :=
  QPochhammer[lambda, q^2] QPochhammer[q^2/lambda, q^2];

(* Keep the original field evaluation and StreamPlot settings. *)
makeDefectPlot[phi_, alpha_ : 0] :=
 Module[{phiPlus, phiMinus, zPlus, zMinus, zetaPlus, zetaMinus,
   twoIOmega, Omega, rawTheta, thetaPlot, nfieldDef}, phiPlus = phi;
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
    Log[(Pcrowdy[lambdaPlus] Pcrowdy[lambdaTildePlus])/(Pcrowdy[
         lambdaMinus] Pcrowdy[lambdaTildeMinus])]];
  Omega[z_] := -I/2 twoIOmega[z];
  rawTheta[x_, y_] := Re[Omega[x + I y]];
  thetaPlot[x_, y_] :=
   rawTheta[x, y] + alpha Log[Sqrt[x^2 + y^2]]/Log[rho];
  nfieldDef[x_, y_] := {Cos[thetaPlot[x, y]], Sin[thetaPlot[x, y]]};
  StreamPlot[nfieldDef[x, y], {x, -H, H}, {y, -H, H},
   RegionFunction ->
    Function[{x, y, vx, vy}, a <= Sqrt[x^2 + y^2] <= b],
   RegionFillingStyle -> None, RegionBoundaryStyle -> None,
   StreamMarkers -> None, StreamColorFunction -> None, Frame -> False,
    Axes -> False, AspectRatio -> 1, PlotRange -> 1.03 {{-H, H}, {-H, H}},
   PlotRangePadding -> None, ImagePadding -> 0, ImageSize -> 220,
   Epilog -> {Black, Thick, Circle[{0, 0}, a], Black, Thick,
     Circle[{0, 0}, b], Blue, Disk[{Re[zPlus], Im[zPlus]}, 0.05], Red,
      Disk[{Re[zMinus], Im[zMinus]}, 0.05]}]];

phiVals = {Pi/12, Pi/3, 2 Pi/3, 11 Pi/12};
alphaVals = {Pi/4, Pi/2, 3 Pi/4};

pltGridAlpha =
 GraphicsGrid[
  Table[makeDefectPlot[#, alpha] & /@ phiVals, {alpha, alphaVals}],
  Spacings -> {10, 10}];

Export["Fig5.pdf", pltGridAlpha, "PDF"];

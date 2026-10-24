(* Figure 11: Defect-mediated unwinding in a confocal elliptic annulus.

   Manuscript coordinates: z = c Exp[I alpha] Cosh[W], W = u + I v.
   The code uses an auxiliary circular coordinate w = Exp[W - uOuter].
   Both boundary ellipses rotate through alpha; director angles are
   measured in the fixed laboratory frame.

   Model: one-constant Frank elasticity with a prescribed +/-1/2 pair.
   The physical anchoring is theta = alpha (mod Pi) on the inclusion
   and theta = 0 (mod Pi) on the outer ellipse.

   Dependencies: Mathematica only; MaTeX is not required.
   Output: Fig11.pdf in Mathematica's current working directory,
   Directory[], which need not be the directory containing this script.
*)

ClearAll[x, y, t, u, w, phi, alpha, lambda,
  c, uInner, uOuter, rho, L, q, H, defectRadius,
  rPlus, rMinus, uPlus, uMinus, phiVals, alphaVals, pltGridAlpha,
  uCoordRot, vCoordRot, vFromLabPhi, wFromXY, zFromW,
  Zmap, Pcrowdy, makeDefectPlot];

(* Geometry. c is the center-to-focus distance; the laboratory foci
   are +/-c Exp[I alpha]. L is the manuscript's ell = uOuter - uInner.
   rho is the INNER RADIUS OF THE AUXILIARY CIRCULAR ANNULUS, not a
   ratio of physical semiaxes. Require 0 < uInner < uOuter and c > 0. *)
c = 1;
uInner = 0.1;
uOuter = 2.0;

rho = Exp[-(uOuter - uInner)];
L = Log[1/rho];
q = Exp[-Pi^2/L];

(* Preserve the original 8% plot margin and graphical marker size.
   defectRadius is a display radius, not the physical core cutoff. *)
H = 1.08 c Cosh[uOuter];
defectRadius = 0.035 H;

(* Fixed defect radii in the AUXILIARY annulus, rho < rPlus,rMinus < 1.
   These are not physical distances from the center. They correspond
   to fixed confocal ellipses uPlus and uMinus, not optimized paths.
   Equal radii give the symmetric pair used in the manuscript. *)
rPlus = 0.7;
rMinus = 0.7;

uPlus = uOuter + Log[rPlus];
uMinus = uOuter + Log[rMinus];

(* Map from auxiliary w to the upper-half-plane coordinate zeta.
   Despite its legacy name, Zmap returns zeta, not the strip coordinate.
   Since Log[w] = W - uOuter on a chosen covering-strip branch,

     Exp[-I Pi Log[w]/L] = -Exp[-I Pi (W - uInner)/L].

   Thus this is exactly the manuscript's map: the minus sign is
   supplied by the shift uOuter - uInner = L. *)
Zmap[w_] := Exp[-I Pi Log[w]/L];

(* Exact prime-function product; no finite-product approximation.
   P(q^2 lambda) = -P(lambda)/lambda ensures periodic director data. *)
Pcrowdy[lambda_] :=
  QPochhammer[lambda, q^2] QPochhammer[q^2/lambda, q^2];

(* Inverse elliptic coordinates after undoing the laboratory rotation.
   The principal v jumps by 2 Pi on the negative major-axis ray, but
   Exp[I v] is continuous there. The focal segment is inside the hole.
   Keep the two coordinate accessors explicit for readability. *)
uCoordRot[x_, y_, alpha_] := Module[{X, Y},
  X = x Cos[alpha] + y Sin[alpha];
  Y = -x Sin[alpha] + y Cos[alpha];
  Re[ArcCosh[(X + I Y)/c]]
];

vCoordRot[x_, y_, alpha_] := Module[{X, Y, xi},
  X = x Cos[alpha] + y Sin[alpha];
  Y = -x Sin[alpha] + y Cos[alpha];
  xi = ArcCosh[(X + I Y)/c];
  Im[xi]
];

(* Convert a LABORATORY polar angle to the elliptic angle v.
   ArcTan[x,y] retains the quadrant; ArcTan[Tan[beta]/Tanh[u]] would
   not. In general v differs from beta = phi - alpha. *)
vFromLabPhi[u_, phi_, alpha_] := Module[{beta},
  beta = phi - alpha;
  ArcTan[Cos[beta]/Cosh[u], Sin[beta]/Sinh[u]]
];

(* Physical position -> auxiliary circular annulus. *)
wFromXY[x_, y_, alpha_] := Module[{u, v},
  u = uCoordRot[x, y, alpha];
  v = vCoordRot[x, y, alpha];
  Exp[(u - uOuter) + I v]
];

(* Auxiliary circular annulus -> physical position. *)
zFromW[w_, alpha_] := Module[{u, v, z0},
  u = uOuter + Log[Abs[w]];
  v = Arg[w];
  z0 = c Cosh[u + I v];
  Exp[I alpha] z0
];

(* A panel at prescribed defect positions.
   phi is the polar half-separation measured from the positive major
   axis, so the laboratory angles are alpha +/- phi. The elliptic
   half-separation is vPlus, not phi. Use 0 < phi < Pi: the coincident
   endpoints require limits that preserve the defects' winding paths. *)
makeDefectPlot[phi_, alpha_ : 0] :=
 Module[{phiPlus, phiMinus, vPlus, vMinus, wPlus, wMinus, zPlus,
   zMinus, zetaPlus, zetaMinus, twoIOmega, Omega, rawThetaD,
   thetaBase, thetaPlot, nfieldDef, p},

  phiPlus = alpha + phi;
  phiMinus = alpha - phi;
  vPlus = vFromLabPhi[uPlus, phiPlus, alpha];
  vMinus = vFromLabPhi[uMinus, phiMinus, alpha];
  wPlus = rPlus Exp[I vPlus];
  wMinus = rMinus Exp[I vMinus];
  zPlus = zFromW[wPlus, alpha];
  zMinus = zFromW[wMinus, alpha];
  zetaPlus = Zmap[wPlus];
  zetaMinus = Zmap[wMinus];

  (* Defect potential and same-sign reflected images. Its real part
     vanishes on both boundaries. Principal-Log jumps change thetaD
     by Pi and therefore leave the physical nematic director unchanged. *)
  twoIOmega[w_] :=
   Module[{zeta, lambdaPlus, lambdaMinus, lambdaTildePlus,
     lambdaTildeMinus},
    zeta = Zmap[w];
    lambdaPlus = zetaPlus/zeta;
    lambdaMinus = zetaMinus/zeta;
    lambdaTildePlus = Conjugate[zetaPlus]/zeta;
    lambdaTildeMinus = Conjugate[zetaMinus]/zeta;
    Log[(Pcrowdy[lambdaPlus] Pcrowdy[lambdaTildePlus])/
      (Pcrowdy[lambdaMinus] Pcrowdy[lambdaTildeMinus])]
   ];
  Omega[w_] := -I/2 twoIOmega[w];
  rawThetaD[x_, y_] := Re[Omega[wFromXY[x, y, alpha]]];

  (* Add the background spiral in the LABORATORY frame, with thetaBase
     equal to alpha on the inclusion and zero on the outer ellipse.
     Do not reduce alpha modulo Pi or add a conformal-map rotation.
     As phi -> Pi from below (equal radii), the total field approaches
     (alpha - Pi) (uOuter - u)/L modulo Pi, consistent with unwinding. *)
  thetaBase[x_, y_] :=
   alpha (uOuter - uCoordRot[x, y, alpha])/(uOuter - uInner);
  thetaPlot[x_, y_] := rawThetaD[x, y] + thetaBase[x, y];
  nfieldDef[x_, y_] := {Cos[thetaPlot[x, y]], Sin[thetaPlot[x, y]]};

  (* Original streamline options, boundary styles, and marker colors. *)
  p = StreamPlot[nfieldDef[x, y], {x, -H, H}, {y, -H, H},
    RegionFunction ->
     Function[{x, y, vx, vy},
      uInner <= uCoordRot[x, y, alpha] <= uOuter],
    RegionFillingStyle -> None, RegionBoundaryStyle -> None,
    StreamMarkers -> None, StreamColorFunction -> None,
    Frame -> False, Axes -> False, AspectRatio -> 1,
    PlotRange -> {{-H, H}, {-H, H}}, PlotRangePadding -> None,
    ImagePadding -> 0, ImageSize -> 220];

  Show[p,
   ParametricPlot[{{Re[Exp[I alpha] c Cosh[uInner + I t]],
      Im[Exp[I alpha] c Cosh[uInner + I t]]}, {Re[
       Exp[I alpha] c Cosh[uOuter + I t]],
      Im[Exp[I alpha] c Cosh[uOuter + I t]]}}, {t, -Pi, Pi},
    PlotStyle -> {Directive[Black, Thick], Directive[Black, Thick]},
    Axes -> False, Frame -> False],
   Graphics[{Blue, Disk[{Re[zPlus], Im[zPlus]}, defectRadius], Red,
     Disk[{Re[zMinus], Im[zMinus]}, defectRadius]}],
   PlotRange -> {{-H, H}, {-H, H}}, Frame -> False, Axes -> False,
   AspectRatio -> 1, ImagePadding -> 5]
 ];

(* Columns: increasing polar half-separation. Rows: inclusion rotation.
   At alpha = Pi/2 the pair starts above the vertical inclusion and
   finishes below it. Horizontal reflection pairs columns 1 with 4
   and 2 with 3, preserving the charges of the reflected defects. *)
phiVals = {Pi/12, Pi/3, 2 Pi/3, 11 Pi/12};
alphaVals = {Pi/4, Pi/2, 3 Pi/4};

pltGridAlpha = GraphicsGrid[
  Table[makeDefectPlot[phi, alpha], {alpha, alphaVals}, {phi, phiVals}],
  Spacings -> {10, 10}
];

Export["Fig11.pdf", pltGridAlpha, "PDF", ImageResolution -> 150];

(*
	LanHEP output produced at Tue Apr 25 17:16:13 2017
	from the file '/export/home/ab1u06/export/belyaev/Dropbox/VDM/vdm_model/vdm_for_fa.lhep'
	Model named 'VDM'
*)


IndexRange[ Index[Colour] ] = NoUnfold[Range[3]]
IndexRange[ Index[Gluon] ] = NoUnfold[Range[8]]

VSESign := -1

		(* Model particles  *)

M$ClassesDescription = {

  U[1] ==  { (* A.c/A.C *)
	SelfConjugate -> False,
	Indices -> {},
	Mass -> 0,
	Charge -> 0,
	PropagatorLabel -> "A.c",
	PropagatorType -> GhostDash,
	PropagatorArrow -> Forward },

  V[1] ==  { (* photon *)
	SelfConjugate -> True,
	Indices -> {},
	Mass -> 0,
	Charge -> 0,
	PropagatorLabel -> "A",
	PropagatorType -> Sine,
	PropagatorArrow -> None },

  U[2] ==  { (* Z.c/Z.C *)
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MZ,
	Charge -> 0,
	PropagatorLabel -> "Z.c",
	PropagatorType -> GhostDash,
	PropagatorArrow -> Forward },

  V[2] ==  { (* 'Z boson' *)
	SelfConjugate -> True,
	Indices -> {},
	Mass -> MZ,
	Charge -> 0,
	PropagatorLabel -> "Z",
	PropagatorType -> Sine,
	PropagatorArrow -> None },

  S[1] ==  { (* 'Z.f' *)
	SelfConjugate -> True,
	Indices -> {},
	Mass -> MZ,
	Charge -> 0,
	PropagatorLabel -> "Z.f",
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None },

  U[3] ==  { (* W+.c/W-.C *)
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MW,
	Charge -> 1,
	PropagatorLabel -> "W+.c",
	PropagatorType -> GhostDash,
	PropagatorArrow -> Forward },

  U[4] ==  { (* W-.c/W+.C *)
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MW,
	Charge -> -1,
	PropagatorLabel -> "W+.c",
	PropagatorType -> GhostDash,
	PropagatorArrow -> Forward },

  V[3] ==  { (* 'W boson' *)
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MW,
	Charge -> 1,
	PropagatorLabel -> "W+",
	PropagatorType -> Sine,
	PropagatorArrow -> Forward },

  S[2] ==  { (* 'W+.f' *)
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MW,
	Charge -> 1,
	PropagatorLabel -> "W+.f",
	PropagatorType -> ScalarDash,
	PropagatorArrow -> Forward },

  U[5] ==  { (* G.c/G.C *)
	SelfConjugate -> False,
	Indices -> {Index[Gluon]},
	Mass -> 0,
	Charge -> 0,
	PropagatorLabel -> "G.c",
	PropagatorType -> GhostDash,
	PropagatorArrow -> Forward },

  V[4] ==  { (* gluon *)
	SelfConjugate -> True,
	Indices -> {Index[Gluon]},
	Mass -> 0,
	Charge -> 0,
	PropagatorLabel -> "G",
	PropagatorType -> Cycles,
	PropagatorArrow -> None },

  V[5] ==  { (* cDM *)
	SelfConjugate -> False,
	Indices -> {},
	Mass -> MVp,
	Charge -> 1,
	PropagatorLabel -> "~V+",
	PropagatorType -> Sine,
	PropagatorArrow -> Forward },

  V[6] ==  { (* nDM *)
	SelfConjugate -> True,
	Indices -> {},
	Mass -> MV0,
	Charge -> 0,
	PropagatorLabel -> "~V0",
	PropagatorType -> Sine,
	PropagatorArrow -> None },

  S[3] ==  { (* Higgs *)
	SelfConjugate -> True,
	Indices -> {},
	Mass -> MH,
	Charge -> 0,
	PropagatorLabel -> "H",
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None },

  F[1] ==  { (* 'u-quark' *)
	SelfConjugate -> False,
	Indices -> {Index[Colour]},
	Mass -> Mu,
	Charge -> 2/3,
	PropagatorLabel -> "u",
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

  F[2] ==  { (* 'd-quark' *)
	SelfConjugate -> False,
	Indices -> {Index[Colour]},
	Mass -> Md,
	Charge -> -1/3,
	PropagatorLabel -> "d",
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

  F[3] ==  { (* 'c-quark' *)
	SelfConjugate -> False,
	Indices -> {Index[Colour]},
	Mass -> Mc,
	Charge -> 2/3,
	PropagatorLabel -> "c",
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

  F[4] ==  { (* 's-quark' *)
	SelfConjugate -> False,
	Indices -> {Index[Colour]},
	Mass -> Ms,
	Charge -> -1/3,
	PropagatorLabel -> "s",
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

  F[5] ==  { (* 't-quark' *)
	SelfConjugate -> False,
	Indices -> {Index[Colour]},
	Mass -> Mtop,
	Charge -> 2/3,
	PropagatorLabel -> "t",
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

  F[6] ==  { (* 'b-quark' *)
	SelfConjugate -> False,
	Indices -> {Index[Colour]},
	Mass -> Mb,
	Charge -> -1/3,
	PropagatorLabel -> "b",
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

  F[7] ==  { (* neutrino *)
	SelfConjugate -> False,
	Indices -> {},
	Mass -> 0,
	Charge -> 0,
	PropagatorLabel -> "ne",
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

  F[8] ==  { (* electron *)
	SelfConjugate -> False,
	Indices -> {},
	Mass -> Me,
	Charge -> -1,
	PropagatorLabel -> "e",
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

  F[9] ==  { (* 'mu-neutrino' *)
	SelfConjugate -> False,
	Indices -> {},
	Mass -> 0,
	Charge -> 0,
	PropagatorLabel -> "nm",
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

  F[10] ==  { (* muon *)
	SelfConjugate -> False,
	Indices -> {},
	Mass -> Mm,
	Charge -> -1,
	PropagatorLabel -> "m",
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

  F[11] ==  { (* 'tau-neutrino' *)
	SelfConjugate -> False,
	Indices -> {},
	Mass -> 0,
	Charge -> 0,
	PropagatorLabel -> "nl",
	PropagatorType -> Straight,
	PropagatorArrow -> Forward },

  F[12] ==  { (* 'tau-lepton' *)
	SelfConjugate -> False,
	Indices -> {},
	Mass -> Mtau,
	Charge -> -1,
	PropagatorLabel -> "l",
	PropagatorType -> Straight,
	PropagatorArrow -> Forward }}

prt["A"] = V[1]
prt["A.c"] = U[1]
prt["A.C"] = -U[1]
prt["Z"] = V[2]
prt["Z.f"] = S[1]
prt["Z.c"] = U[2]
prt["Z.C"] = -U[2]
prt["W-"] = -V[3]
prt["W-.f"] = -S[2]
prt["W-.c"] = U[4]
prt["W-.C"] = -U[3]
prt["W+"] = V[3]
prt["W+.f"] = S[2]
prt["W+.c"] = U[3]
prt["W+.C"] = -U[4]
prt["G"] = V[4]
prt["G.c"] = U[5]
prt["G.C"] = -U[5]
prt["~V-"] = -V[5]
prt["~V+"] = V[5]
prt["~V0"] = V[6]
prt["H"] = S[3]
prt["U"] = -F[1]
prt["u"] = F[1]
prt["D"] = -F[2]
prt["d"] = F[2]
prt["C"] = -F[3]
prt["c"] = F[3]
prt["S"] = -F[4]
prt["s"] = F[4]
prt["T"] = -F[5]
prt["t"] = F[5]
prt["B"] = -F[6]
prt["b"] = F[6]
prt["Ne"] = -F[7]
prt["ne"] = F[7]
prt["E"] = -F[8]
prt["e"] = F[8]
prt["Nm"] = -F[9]
prt["nm"] = F[9]
prt["M"] = -F[10]
prt["m"] = F[10]
prt["Nl"] = -F[11]
prt["nl"] = F[11]
prt["L"] = -F[12]
prt["l"] = F[12]

Mu[_] = Mu
Md[_] = Md
Mc[_] = Mc
Ms[_] = Ms
Mtop[_] = Mtop
Mb[_] = Mb


FAGaugeXi[_] = 1


M$CouplingMatrices = {

  (*------    W+.C  A.c  W-.f  ------*)
   C[ -U[4], U[1], -S[2] ] == EE MW *
{ 
 { 1 } 
},
  (*------    W+.C  W-.c  H  ------*)
   C[ -U[4], U[4], S[3] ] == -1/2 I EE MW / SW *
{ 
 { 1 } 
},
  (*------    W+.C  W-.c  Z.f  ------*)
   C[ -U[4], U[4], S[1] ] == -1/2 EE MW / SW *
{ 
 { 1 } 
},
  (*------    W+.C  Z.c  W-.f  ------*)
   C[ -U[4], U[2], -S[2] ] == 1/2 EE MW / SW *
{ 
 { 2 CW -1 / CW } 
},
  (*------    W-.C  A.c  W+.f  ------*)
   C[ -U[3], U[1], S[2] ] == - EE MW *
{ 
 { 1 } 
},
  (*------    W-.C  W+.c  H  ------*)
   C[ -U[3], U[3], S[3] ] == -1/2 I EE MW / SW *
{ 
 { 1 } 
},
  (*------    W-.C  W+.c  Z.f  ------*)
   C[ -U[3], U[3], S[1] ] == 1/2 EE MW / SW *
{ 
 { 1 } 
},
  (*------    W-.C  Z.c  W+.f  ------*)
   C[ -U[3], U[2], S[2] ] == -1/2 EE MW / SW *
{ 
 { 2 CW -1 / CW } 
},
  (*------    Z.C  W+.c  W-.f  ------*)
   C[ -U[2], U[3], -S[2] ] == -1/2 / CW EE MW / SW *
{ 
 { 1 } 
},
  (*------    Z.C  W-.c  W+.f  ------*)
   C[ -U[2], U[4], S[2] ] == 1/2 / CW EE MW / SW *
{ 
 { 1 } 
},
  (*------    Z.C  Z.c  H  ------*)
   C[ -U[2], U[2], S[3] ] == -1/2 I / CW^2 EE MW / SW *
{ 
 { 1 } 
},
  (*------    A.C  W+.c  W-  ------*)
   C[ -U[1], U[3], -V[3] ] == -2 I EE *
{ 
 { 1 },
 { 0 } 
},
  (*------    A.C  W-.c  W+  ------*)
   C[ -U[1], U[4], V[3] ] == 2 I EE *
{ 
 { 1 },
 { 0 } 
},
  (*------    G.C  G.c  G  ------*)
   C[ -U[5,{c1}], U[5,{c2}], V[4,{c3}] ] == - GG FASUNF[c1, c2, c3] *
{ 
 { 0 },
 { 1 } 
},
  (*------    W+.C  A.c  W-  ------*)
   C[ -U[4], U[1], -V[3] ] == 2 I EE *
{ 
 { 1 },
 { 0 } 
},
  (*------    W+.C  W-.c  A  ------*)
   C[ -U[4], U[4], V[1] ] == -2 I EE *
{ 
 { 1 },
 { 0 } 
},
  (*------    W+.C  W-.c  Z  ------*)
   C[ -U[4], U[4], V[2] ] == -2 I CW EE / SW *
{ 
 { 1 },
 { 0 } 
},
  (*------    W+.C  Z.c  W-  ------*)
   C[ -U[4], U[2], -V[3] ] == 2 I CW EE / SW *
{ 
 { 1 },
 { 0 } 
},
  (*------    W-.C  A.c  W+  ------*)
   C[ -U[3], U[1], V[3] ] == -2 I EE *
{ 
 { 1 },
 { 0 } 
},
  (*------    W-.C  W+.c  A  ------*)
   C[ -U[3], U[3], V[1] ] == 2 I EE *
{ 
 { 1 },
 { 0 } 
},
  (*------    W-.C  W+.c  Z  ------*)
   C[ -U[3], U[3], V[2] ] == 2 I CW EE / SW *
{ 
 { 1 },
 { 0 } 
},
  (*------    W-.C  Z.c  W+  ------*)
   C[ -U[3], U[2], V[3] ] == -2 I CW EE / SW *
{ 
 { 1 },
 { 0 } 
},
  (*------    Z.C  W+.c  W-  ------*)
   C[ -U[2], U[3], -V[3] ] == -2 I CW EE / SW *
{ 
 { 1 },
 { 0 } 
},
  (*------    Z.C  W-.c  W+  ------*)
   C[ -U[2], U[4], V[3] ] == 2 I CW EE / SW *
{ 
 { 1 },
 { 0 } 
},
  (*------    H  H  H  ------*)
   C[ S[3], S[3], S[3] ] == -12 I / EE MW SW lam *
{ 
 { 1 } 
},
  (*------    H  W+.f  W-.f  ------*)
   C[ S[3], S[2], -S[2] ] == -4 I / EE MW SW lam *
{ 
 { 1 } 
},
  (*------    H  Z.f  Z.f  ------*)
   C[ S[3], S[1], S[1] ] == -4 I / EE MW SW lam *
{ 
 { 1 } 
},
  (*------    H  W+.f  W-  ------*)
   C[ S[3], S[2], -V[3] ] == 1/2 EE / SW *
{ 
 { 1 },
 { -1 } 
},
  (*------    H  W-.f  W+  ------*)
   C[ S[3], -S[2], V[3] ] == -1/2 EE / SW *
{ 
 { -1 },
 { 1 } 
},
  (*------    H  Z.f  Z  ------*)
   C[ S[3], S[1], V[2] ] == 1/2 / CW EE / SW *
{ 
 { 1 },
 { -1 } 
},
  (*------    W+.f  W-.f  A  ------*)
   C[ S[2], -S[2], V[1] ] == I EE *
{ 
 { -1 },
 { 1 } 
},
  (*------    W+.f  W-.f  Z  ------*)
   C[ S[2], -S[2], V[2] ] == 1/2 I EE / SW *
{ 
 { -2 CW + 1 / CW },
 { 2 CW -1 / CW } 
},
  (*------    W+.f  Z.f  W-  ------*)
   C[ S[2], S[1], -V[3] ] == -1/2 I EE / SW *
{ 
 { -1 },
 { 1 } 
},
  (*------    W-.f  Z.f  W+  ------*)
   C[ -S[2], S[1], V[3] ] == -1/2 I EE / SW *
{ 
 { 1 },
 { -1 } 
},
  (*------    H  W+  W-  ------*)
   C[ S[3], V[3], -V[3] ] == I EE MW / SW *
{ 
 { 1 } 
},
  (*------    H  Z  Z  ------*)
   C[ S[3], V[2], V[2] ] == I / CW^2 EE MW / SW *
{ 
 { 1 } 
},
  (*------    H  ~V+  ~V-  ------*)
   C[ S[3], V[5], -V[5] ] == 4 I / EE MW SW a *
{ 
 { 1 } 
},
  (*------    H  ~V0  ~V0  ------*)
   C[ S[3], V[6], V[6] ] == 4 I / EE MW SW a *
{ 
 { 1 } 
},
  (*------    W+.f  A  W-  ------*)
   C[ S[2], V[1], -V[3] ] == EE MW *
{ 
 { 1 } 
},
  (*------    W+.f  W-  Z  ------*)
   C[ S[2], -V[3], V[2] ] == -1 / CW EE MW SW *
{ 
 { 1 } 
},
  (*------    W-.f  A  W+  ------*)
   C[ -S[2], V[1], V[3] ] == - EE MW *
{ 
 { 1 } 
},
  (*------    W-.f  W+  Z  ------*)
   C[ -S[2], V[3], V[2] ] == 1 / CW EE MW SW *
{ 
 { 1 } 
},
  (*------    B  b  H  ------*)
   C[ -F[6,{c1}], F[6,{c2}], S[3] ] == -1/2 I EE / MW Mb / SW IndexDelta[c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    B  b  Z.f  ------*)
   C[ -F[6,{c1}], F[6,{c2}], S[1] ] == 1/2 EE / MW Mb / SW IndexDelta[c1, c2] *
{ 
 { 1 },
 { -1 } 
},
  (*------    B  c  W-.f  ------*)
   C[ -F[6,{c1}], F[3,{c2}], -S[2] ] == 1/2 EE / MW / SW Sqrt2 Vcb *
		IndexDelta[c1, c2] *
{ 
 { - Mc },
 { Mb } 
},
  (*------    B  t  W-.f  ------*)
   C[ -F[6,{c1}], F[5,{c2}], -S[2] ] == 1/2 EE / MW / SW Sqrt2 Vtb *
		IndexDelta[c1, c2] *
{ 
 { - Mtop },
 { Mb } 
},
  (*------    B  u  W-.f  ------*)
   C[ -F[6,{c1}], F[1,{c2}], -S[2] ] == 1/2 EE / MW / SW Sqrt2 Vub *
		IndexDelta[c1, c2] *
{ 
 { - Mu },
 { Mb } 
},
  (*------    C  b  W+.f  ------*)
   C[ -F[3,{c1}], F[6,{c2}], S[2] ] == -1/2 EE / MW / SW Sqrt2 Vcb *
		IndexDelta[c1, c2] *
{ 
 { Mb },
 { - Mc } 
},
  (*------    C  c  H  ------*)
   C[ -F[3,{c1}], F[3,{c2}], S[3] ] == -1/2 I EE / MW Mc / SW IndexDelta[c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    C  c  Z.f  ------*)
   C[ -F[3,{c1}], F[3,{c2}], S[1] ] == -1/2 EE / MW Mc / SW IndexDelta[c1, c2] *
{ 
 { 1 },
 { -1 } 
},
  (*------    C  d  W+.f  ------*)
   C[ -F[3,{c1}], F[2,{c2}], S[2] ] == -1/2 EE / MW / SW Sqrt2 Vcd *
		IndexDelta[c1, c2] *
{ 
 { Md },
 { - Mc } 
},
  (*------    C  s  W+.f  ------*)
   C[ -F[3,{c1}], F[4,{c2}], S[2] ] == -1/2 EE / MW / SW Sqrt2 Vcs *
		IndexDelta[c1, c2] *
{ 
 { Ms },
 { - Mc } 
},
  (*------    D  c  W-.f  ------*)
   C[ -F[2,{c1}], F[3,{c2}], -S[2] ] == 1/2 EE / MW / SW Sqrt2 Vcd *
		IndexDelta[c1, c2] *
{ 
 { - Mc },
 { Md } 
},
  (*------    D  d  H  ------*)
   C[ -F[2,{c1}], F[2,{c2}], S[3] ] == -1/2 I EE / MW Md / SW IndexDelta[c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    D  d  Z.f  ------*)
   C[ -F[2,{c1}], F[2,{c2}], S[1] ] == 1/2 EE / MW Md / SW IndexDelta[c1, c2] *
{ 
 { 1 },
 { -1 } 
},
  (*------    D  t  W-.f  ------*)
   C[ -F[2,{c1}], F[5,{c2}], -S[2] ] == 1/2 EE / MW / SW Sqrt2 Vtd *
		IndexDelta[c1, c2] *
{ 
 { - Mtop },
 { Md } 
},
  (*------    D  u  W-.f  ------*)
   C[ -F[2,{c1}], F[1,{c2}], -S[2] ] == 1/2 EE / MW / SW Sqrt2 Vud *
		IndexDelta[c1, c2] *
{ 
 { - Mu },
 { Md } 
},
  (*------    E  e  H  ------*)
   C[ -F[8], F[8], S[3] ] == -1/2 I EE / MW Me / SW *
{ 
 { 1 },
 { 1 } 
},
  (*------    E  e  Z.f  ------*)
   C[ -F[8], F[8], S[1] ] == 1/2 EE / MW Me / SW *
{ 
 { 1 },
 { -1 } 
},
  (*------    E  ne  W-.f  ------*)
   C[ -F[8], F[7], -S[2] ] == 1/2 EE / MW Me / SW Sqrt2 *
{ 
 { 0 },
 { 1 } 
},
  (*------    L  l  H  ------*)
   C[ -F[12], F[12], S[3] ] == -1/2 I EE / MW Mtau / SW *
{ 
 { 1 },
 { 1 } 
},
  (*------    L  l  Z.f  ------*)
   C[ -F[12], F[12], S[1] ] == 1/2 EE / MW Mtau / SW *
{ 
 { 1 },
 { -1 } 
},
  (*------    L  nl  W-.f  ------*)
   C[ -F[12], F[11], -S[2] ] == 1/2 EE / MW Mtau / SW Sqrt2 *
{ 
 { 0 },
 { 1 } 
},
  (*------    M  m  H  ------*)
   C[ -F[10], F[10], S[3] ] == -1/2 I EE / MW Mm / SW *
{ 
 { 1 },
 { 1 } 
},
  (*------    M  m  Z.f  ------*)
   C[ -F[10], F[10], S[1] ] == 1/2 EE / MW Mm / SW *
{ 
 { 1 },
 { -1 } 
},
  (*------    M  nm  W-.f  ------*)
   C[ -F[10], F[9], -S[2] ] == 1/2 EE / MW Mm / SW Sqrt2 *
{ 
 { 0 },
 { 1 } 
},
  (*------    Ne  e  W+.f  ------*)
   C[ -F[7], F[8], S[2] ] == -1/2 EE / MW Me / SW Sqrt2 *
{ 
 { 1 },
 { 0 } 
},
  (*------    Nl  l  W+.f  ------*)
   C[ -F[11], F[12], S[2] ] == -1/2 EE / MW Mtau / SW Sqrt2 *
{ 
 { 1 },
 { 0 } 
},
  (*------    Nm  m  W+.f  ------*)
   C[ -F[9], F[10], S[2] ] == -1/2 EE / MW Mm / SW Sqrt2 *
{ 
 { 1 },
 { 0 } 
},
  (*------    S  c  W-.f  ------*)
   C[ -F[4,{c1}], F[3,{c2}], -S[2] ] == 1/2 EE / MW / SW Sqrt2 Vcs *
		IndexDelta[c1, c2] *
{ 
 { - Mc },
 { Ms } 
},
  (*------    S  s  H  ------*)
   C[ -F[4,{c1}], F[4,{c2}], S[3] ] == -1/2 I EE / MW Ms / SW IndexDelta[c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    S  s  Z.f  ------*)
   C[ -F[4,{c1}], F[4,{c2}], S[1] ] == 1/2 EE / MW Ms / SW IndexDelta[c1, c2] *
{ 
 { 1 },
 { -1 } 
},
  (*------    S  t  W-.f  ------*)
   C[ -F[4,{c1}], F[5,{c2}], -S[2] ] == 1/2 EE / MW / SW Sqrt2 Vts *
		IndexDelta[c1, c2] *
{ 
 { - Mtop },
 { Ms } 
},
  (*------    S  u  W-.f  ------*)
   C[ -F[4,{c1}], F[1,{c2}], -S[2] ] == 1/2 EE / MW / SW Sqrt2 Vus *
		IndexDelta[c1, c2] *
{ 
 { - Mu },
 { Ms } 
},
  (*------    T  b  W+.f  ------*)
   C[ -F[5,{c1}], F[6,{c2}], S[2] ] == -1/2 EE / MW / SW Sqrt2 Vtb *
		IndexDelta[c1, c2] *
{ 
 { Mb },
 { - Mtop } 
},
  (*------    T  d  W+.f  ------*)
   C[ -F[5,{c1}], F[2,{c2}], S[2] ] == -1/2 EE / MW / SW Sqrt2 Vtd *
		IndexDelta[c1, c2] *
{ 
 { Md },
 { - Mtop } 
},
  (*------    T  s  W+.f  ------*)
   C[ -F[5,{c1}], F[4,{c2}], S[2] ] == -1/2 EE / MW / SW Sqrt2 Vts *
		IndexDelta[c1, c2] *
{ 
 { Ms },
 { - Mtop } 
},
  (*------    T  t  H  ------*)
   C[ -F[5,{c1}], F[5,{c2}], S[3] ] == -1/2 I EE / MW Mtop / SW *
		IndexDelta[c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    T  t  Z.f  ------*)
   C[ -F[5,{c1}], F[5,{c2}], S[1] ] == -1/2 EE / MW Mtop / SW IndexDelta[c1, c2] *
{ 
 { 1 },
 { -1 } 
},
  (*------    U  b  W+.f  ------*)
   C[ -F[1,{c1}], F[6,{c2}], S[2] ] == -1/2 EE / MW / SW Sqrt2 Vub *
		IndexDelta[c1, c2] *
{ 
 { Mb },
 { - Mu } 
},
  (*------    U  d  W+.f  ------*)
   C[ -F[1,{c1}], F[2,{c2}], S[2] ] == -1/2 EE / MW / SW Sqrt2 Vud *
		IndexDelta[c1, c2] *
{ 
 { Md },
 { - Mu } 
},
  (*------    U  s  W+.f  ------*)
   C[ -F[1,{c1}], F[4,{c2}], S[2] ] == -1/2 EE / MW / SW Sqrt2 Vus *
		IndexDelta[c1, c2] *
{ 
 { Ms },
 { - Mu } 
},
  (*------    U  u  H  ------*)
   C[ -F[1,{c1}], F[1,{c2}], S[3] ] == -1/2 I EE / MW Mu / SW IndexDelta[c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    U  u  Z.f  ------*)
   C[ -F[1,{c1}], F[1,{c2}], S[1] ] == -1/2 EE / MW Mu / SW IndexDelta[c1, c2] *
{ 
 { 1 },
 { -1 } 
},
  (*------    B  b  A  ------*)
   C[ -F[6,{c1}], F[6,{c2}], V[1] ] == 1/3 I EE IndexDelta[c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    B  b  G  ------*)
   C[ -F[6,{c1}], F[6,{c2}], V[4,{c3}] ] == - I GG FASUNT[c3, c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    B  b  Z  ------*)
   C[ -F[6,{c1}], F[6,{c2}], V[2] ] == 1/6 I EE IndexDelta[c1, c2] *
{ 
 { 2 CW / SW + 1 / CW / SW },
 { -2 / CW SW } 
},
  (*------    B  t  W-  ------*)
   C[ -F[6,{c1}], F[5,{c2}], -V[3] ] == -1/2 I EE / SW Sqrt2 IndexDelta[c1, c2] *
{ 
 { 1 },
 { 0 } 
},
  (*------    C  c  A  ------*)
   C[ -F[3,{c1}], F[3,{c2}], V[1] ] == -2/3 I EE IndexDelta[c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    C  c  G  ------*)
   C[ -F[3,{c1}], F[3,{c2}], V[4,{c3}] ] == - I GG FASUNT[c3, c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    C  c  Z  ------*)
   C[ -F[3,{c1}], F[3,{c2}], V[2] ] == -1/6 I EE IndexDelta[c1, c2] *
{ 
 { 4 CW / SW -1 / CW / SW },
 { -4 / CW SW } 
},
  (*------    C  s  W+  ------*)
   C[ -F[3,{c1}], F[4,{c2}], V[3] ] == -1/2 I EE / SW Sqrt2 IndexDelta[c1, c2] *
{ 
 { 1 },
 { 0 } 
},
  (*------    D  d  A  ------*)
   C[ -F[2,{c1}], F[2,{c2}], V[1] ] == 1/3 I EE IndexDelta[c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    D  d  G  ------*)
   C[ -F[2,{c1}], F[2,{c2}], V[4,{c3}] ] == - I GG FASUNT[c3, c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    D  d  Z  ------*)
   C[ -F[2,{c1}], F[2,{c2}], V[2] ] == 1/6 I EE IndexDelta[c1, c2] *
{ 
 { 2 CW / SW + 1 / CW / SW },
 { -2 / CW SW } 
},
  (*------    D  u  W-  ------*)
   C[ -F[2,{c1}], F[1,{c2}], -V[3] ] == -1/2 I EE / SW Sqrt2 IndexDelta[c1, c2] *
{ 
 { 1 },
 { 0 } 
},
  (*------    E  e  A  ------*)
   C[ -F[8], F[8], V[1] ] == I EE *
{ 
 { 1 },
 { 1 } 
},
  (*------    E  e  Z  ------*)
   C[ -F[8], F[8], V[2] ] == 1/2 I EE *
{ 
 { 2 CW / SW -1 / CW / SW },
 { -2 / CW SW } 
},
  (*------    E  ne  W-  ------*)
   C[ -F[8], F[7], -V[3] ] == -1/2 I EE / SW Sqrt2 *
{ 
 { 1 },
 { 0 } 
},
  (*------    L  l  A  ------*)
   C[ -F[12], F[12], V[1] ] == I EE *
{ 
 { 1 },
 { 1 } 
},
  (*------    L  l  Z  ------*)
   C[ -F[12], F[12], V[2] ] == 1/2 I EE *
{ 
 { 2 CW / SW -1 / CW / SW },
 { -2 / CW SW } 
},
  (*------    L  nl  W-  ------*)
   C[ -F[12], F[11], -V[3] ] == -1/2 I EE / SW Sqrt2 *
{ 
 { 1 },
 { 0 } 
},
  (*------    M  m  A  ------*)
   C[ -F[10], F[10], V[1] ] == I EE *
{ 
 { 1 },
 { 1 } 
},
  (*------    M  m  Z  ------*)
   C[ -F[10], F[10], V[2] ] == 1/2 I EE *
{ 
 { 2 CW / SW -1 / CW / SW },
 { -2 / CW SW } 
},
  (*------    M  nm  W-  ------*)
   C[ -F[10], F[9], -V[3] ] == -1/2 I EE / SW Sqrt2 *
{ 
 { 1 },
 { 0 } 
},
  (*------    Ne  e  W+  ------*)
   C[ -F[7], F[8], V[3] ] == -1/2 I EE / SW Sqrt2 *
{ 
 { 1 },
 { 0 } 
},
  (*------    Ne  ne  Z  ------*)
   C[ -F[7], F[7], V[2] ] == -1/2 I / CW EE / SW *
{ 
 { 1 },
 { 0 } 
},
  (*------    Nl  l  W+  ------*)
   C[ -F[11], F[12], V[3] ] == -1/2 I EE / SW Sqrt2 *
{ 
 { 1 },
 { 0 } 
},
  (*------    Nl  nl  Z  ------*)
   C[ -F[11], F[11], V[2] ] == -1/2 I / CW EE / SW *
{ 
 { 1 },
 { 0 } 
},
  (*------    Nm  m  W+  ------*)
   C[ -F[9], F[10], V[3] ] == -1/2 I EE / SW Sqrt2 *
{ 
 { 1 },
 { 0 } 
},
  (*------    Nm  nm  Z  ------*)
   C[ -F[9], F[9], V[2] ] == -1/2 I / CW EE / SW *
{ 
 { 1 },
 { 0 } 
},
  (*------    S  c  W-  ------*)
   C[ -F[4,{c1}], F[3,{c2}], -V[3] ] == -1/2 I EE / SW Sqrt2 IndexDelta[c1, c2] *
{ 
 { 1 },
 { 0 } 
},
  (*------    S  s  A  ------*)
   C[ -F[4,{c1}], F[4,{c2}], V[1] ] == 1/3 I EE IndexDelta[c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    S  s  G  ------*)
   C[ -F[4,{c1}], F[4,{c2}], V[4,{c3}] ] == - I GG FASUNT[c3, c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    S  s  Z  ------*)
   C[ -F[4,{c1}], F[4,{c2}], V[2] ] == 1/6 I EE IndexDelta[c1, c2] *
{ 
 { 2 CW / SW + 1 / CW / SW },
 { -2 / CW SW } 
},
  (*------    T  b  W+  ------*)
   C[ -F[5,{c1}], F[6,{c2}], V[3] ] == -1/2 I EE / SW Sqrt2 IndexDelta[c1, c2] *
{ 
 { 1 },
 { 0 } 
},
  (*------    T  t  A  ------*)
   C[ -F[5,{c1}], F[5,{c2}], V[1] ] == -2/3 I EE IndexDelta[c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    T  t  G  ------*)
   C[ -F[5,{c1}], F[5,{c2}], V[4,{c3}] ] == - I GG FASUNT[c3, c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    T  t  Z  ------*)
   C[ -F[5,{c1}], F[5,{c2}], V[2] ] == -1/6 I EE IndexDelta[c1, c2] *
{ 
 { 4 CW / SW -1 / CW / SW },
 { -4 / CW SW } 
},
  (*------    U  d  W+  ------*)
   C[ -F[1,{c1}], F[2,{c2}], V[3] ] == -1/2 I EE / SW Sqrt2 IndexDelta[c1, c2] *
{ 
 { 1 },
 { 0 } 
},
  (*------    U  u  A  ------*)
   C[ -F[1,{c1}], F[1,{c2}], V[1] ] == -2/3 I EE IndexDelta[c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    U  u  G  ------*)
   C[ -F[1,{c1}], F[1,{c2}], V[4,{c3}] ] == - I GG FASUNT[c3, c1, c2] *
{ 
 { 1 },
 { 1 } 
},
  (*------    U  u  Z  ------*)
   C[ -F[1,{c1}], F[1,{c2}], V[2] ] == -1/6 I EE IndexDelta[c1, c2] *
{ 
 { 4 CW / SW -1 / CW / SW },
 { -4 / CW SW } 
},
  (*------    A  W+  W-  ------*)
   C[ V[1], V[3], -V[3] ] == - I EE *
{ 
 { 1 },
 { -1 },
 { -1 },
 { 1 },
 { 1 },
 { -1 } 
},
  (*------    A  ~V+  ~V-  ------*)
   C[ V[1], V[5], -V[5] ] == - I EE *
{ 
 { 1 },
 { -1 },
 { -1 },
 { 1 },
 { 1 },
 { -1 } 
},
  (*------    G  G  G  ------*)
   C[ V[4,{c1}], V[4,{c2}], V[4,{c3}] ] == - GG FASUNF[c1, c2, c3] *
{ 
 { -1 },
 { 1 },
 { 1 },
 { -1 },
 { -1 },
 { 1 } 
},
  (*------    W+  W-  Z  ------*)
   C[ V[3], -V[3], V[2] ] == - I CW EE / SW *
{ 
 { 1 },
 { -1 },
 { -1 },
 { 1 },
 { 1 },
 { -1 } 
},
  (*------    W+  ~V-  ~V0  ------*)
   C[ V[3], -V[5], V[6] ] == - I EE / SW *
{ 
 { 1 },
 { -1 },
 { -1 },
 { 1 },
 { 1 },
 { -1 } 
},
  (*------    W-  ~V+  ~V0  ------*)
   C[ -V[3], V[5], V[6] ] == I EE / SW *
{ 
 { 1 },
 { -1 },
 { -1 },
 { 1 },
 { 1 },
 { -1 } 
},
  (*------    Z  ~V+  ~V-  ------*)
   C[ V[2], V[5], -V[5] ] == - I CW EE / SW *
{ 
 { 1 },
 { -1 },
 { -1 },
 { 1 },
 { 1 },
 { -1 } 
},
  (*------    H  H  H  H  ------*)
   C[ S[3], S[3], S[3], S[3] ] == -6 I lam *
{ 
 { 1 } 
},
  (*------    H  H  W+.f  W-.f  ------*)
   C[ S[3], S[3], S[2], -S[2] ] == -2 I lam *
{ 
 { 1 } 
},
  (*------    H  H  Z.f  Z.f  ------*)
   C[ S[3], S[3], S[1], S[1] ] == -2 I lam *
{ 
 { 1 } 
},
  (*------    W+.f  W+.f  W-.f  W-.f  ------*)
   C[ S[2], S[2], -S[2], -S[2] ] == -4 I lam *
{ 
 { 1 } 
},
  (*------    W+.f  W-.f  Z.f  Z.f  ------*)
   C[ S[2], -S[2], S[1], S[1] ] == -2 I lam *
{ 
 { 1 } 
},
  (*------    Z.f  Z.f  Z.f  Z.f  ------*)
   C[ S[1], S[1], S[1], S[1] ] == -6 I lam *
{ 
 { 1 } 
},
  (*------    H  H  W+  W-  ------*)
   C[ S[3], S[3], V[3], -V[3] ] == 1/2 I EE^2 / SW^2 *
{ 
 { 1 } 
},
  (*------    H  H  Z  Z  ------*)
   C[ S[3], S[3], V[2], V[2] ] == 1/2 I / CW^2 EE^2 / SW^2 *
{ 
 { 1 } 
},
  (*------    H  H  ~V+  ~V-  ------*)
   C[ S[3], S[3], V[5], -V[5] ] == 2 I a *
{ 
 { 1 } 
},
  (*------    H  H  ~V0  ~V0  ------*)
   C[ S[3], S[3], V[6], V[6] ] == 2 I a *
{ 
 { 1 } 
},
  (*------    H  W+.f  A  W-  ------*)
   C[ S[3], S[2], V[1], -V[3] ] == 1/2 EE^2 / SW *
{ 
 { 1 } 
},
  (*------    H  W+.f  W-  Z  ------*)
   C[ S[3], S[2], -V[3], V[2] ] == -1/2 / CW EE^2 *
{ 
 { 1 } 
},
  (*------    H  W-.f  A  W+  ------*)
   C[ S[3], -S[2], V[1], V[3] ] == -1/2 EE^2 / SW *
{ 
 { 1 } 
},
  (*------    H  W-.f  W+  Z  ------*)
   C[ S[3], -S[2], V[3], V[2] ] == 1/2 / CW EE^2 *
{ 
 { 1 } 
},
  (*------    W+.f  W-.f  A  A  ------*)
   C[ S[2], -S[2], V[1], V[1] ] == 2 I EE^2 *
{ 
 { 1 } 
},
  (*------    W+.f  W-.f  A  Z  ------*)
   C[ S[2], -S[2], V[1], V[2] ] == I EE^2 / SW *
{ 
 { 2 CW -1 / CW } 
},
  (*------    W+.f  W-.f  W+  W-  ------*)
   C[ S[2], -S[2], V[3], -V[3] ] == 1/2 I EE^2 / SW^2 *
{ 
 { 1 } 
},
  (*------    W+.f  W-.f  Z  Z  ------*)
   C[ S[2], -S[2], V[2], V[2] ] == 1/2 I EE^2 / SW^2 *
{ 
 { 1 / CW^2 + 4 CW^2 -4 } 
},
  (*------    W+.f  W-.f  ~V+  ~V-  ------*)
   C[ S[2], -S[2], V[5], -V[5] ] == 2 I a *
{ 
 { 1 } 
},
  (*------    W+.f  W-.f  ~V0  ~V0  ------*)
   C[ S[2], -S[2], V[6], V[6] ] == 2 I a *
{ 
 { 1 } 
},
  (*------    W+.f  Z.f  A  W-  ------*)
   C[ S[2], S[1], V[1], -V[3] ] == -1/2 I EE^2 / SW *
{ 
 { 1 } 
},
  (*------    W+.f  Z.f  W-  Z  ------*)
   C[ S[2], S[1], -V[3], V[2] ] == 1/2 I / CW EE^2 *
{ 
 { 1 } 
},
  (*------    W-.f  Z.f  A  W+  ------*)
   C[ -S[2], S[1], V[1], V[3] ] == -1/2 I EE^2 / SW *
{ 
 { 1 } 
},
  (*------    W-.f  Z.f  W+  Z  ------*)
   C[ -S[2], S[1], V[3], V[2] ] == 1/2 I / CW EE^2 *
{ 
 { 1 } 
},
  (*------    Z.f  Z.f  W+  W-  ------*)
   C[ S[1], S[1], V[3], -V[3] ] == 1/2 I EE^2 / SW^2 *
{ 
 { 1 } 
},
  (*------    Z.f  Z.f  Z  Z  ------*)
   C[ S[1], S[1], V[2], V[2] ] == 1/2 I / CW^2 EE^2 / SW^2 *
{ 
 { 1 } 
},
  (*------    Z.f  Z.f  ~V+  ~V-  ------*)
   C[ S[1], S[1], V[5], -V[5] ] == 2 I a *
{ 
 { 1 } 
},
  (*------    Z.f  Z.f  ~V0  ~V0  ------*)
   C[ S[1], S[1], V[6], V[6] ] == 2 I a *
{ 
 { 1 } 
},
  (*------    A  A  W+  W-  ------*)
   C[ V[1], V[1], V[3], -V[3] ] == - I EE^2 *
{ 
 { 2 },
 { -1 },
 { -1 } 
},
  (*------    A  A  ~V+  ~V-  ------*)
   C[ V[1], V[1], V[5], -V[5] ] == - I EE^2 *
{ 
 { 2 },
 { -1 },
 { -1 } 
},
  (*------    A  W+  W-  Z  ------*)
   C[ V[1], V[3], -V[3], V[2] ] == - I CW EE^2 / SW *
{ 
 { -1 },
 { -1 },
 { 2 } 
},
  (*------    A  W+  ~V-  ~V0  ------*)
   C[ V[1], V[3], -V[5], V[6] ] == - I EE^2 / SW *
{ 
 { -1 },
 { -1 },
 { 2 } 
},
  (*------    A  W-  ~V+  ~V0  ------*)
   C[ V[1], -V[3], V[5], V[6] ] == - I EE^2 / SW *
{ 
 { -1 },
 { -1 },
 { 2 } 
},
  (*------    A  Z  ~V+  ~V-  ------*)
   C[ V[1], V[2], V[5], -V[5] ] == - I CW EE^2 / SW *
{ 
 { 2 },
 { -1 },
 { -1 } 
},
  (*------    G  G  G  G  ------*)
   C[ V[4,{c1}], V[4,{c2}], V[4,{c3}], V[4,{c4}] ] == - I GG^2 *
{ 
 { FASUNF[c1, c3, c2, c4] + FASUNF[c1, c4, c2, c3] },
 { FASUNF[c1, c2, c3, c4] - FASUNF[c1, c4, c2, c3] },
 { - FASUNF[c1, c2, c3, c4] - FASUNF[c1, c3, c2, c4] } 
},
  (*------    W+  W+  W-  W-  ------*)
   C[ V[3], V[3], -V[3], -V[3] ] == I EE^2 / SW^2 *
{ 
 { 2 },
 { -1 },
 { -1 } 
},
  (*------    W+  W+  ~V-  ~V-  ------*)
   C[ V[3], V[3], -V[5], -V[5] ] == I EE^2 / SW^2 *
{ 
 { 2 },
 { -1 },
 { -1 } 
},
  (*------    W+  W-  Z  Z  ------*)
   C[ V[3], -V[3], V[2], V[2] ] == - I CW^2 EE^2 / SW^2 *
{ 
 { 2 },
 { -1 },
 { -1 } 
},
  (*------    W+  W-  ~V+  ~V-  ------*)
   C[ V[3], -V[3], V[5], -V[5] ] == I EE^2 / SW^2 *
{ 
 { -1 },
 { 2 },
 { -1 } 
},
  (*------    W+  W-  ~V0  ~V0  ------*)
   C[ V[3], -V[3], V[6], V[6] ] == - I EE^2 / SW^2 *
{ 
 { 2 },
 { -1 },
 { -1 } 
},
  (*------    W+  Z  ~V-  ~V0  ------*)
   C[ V[3], V[2], -V[5], V[6] ] == - I CW EE^2 / SW^2 *
{ 
 { -1 },
 { 2 },
 { -1 } 
},
  (*------    W-  W-  ~V+  ~V+  ------*)
   C[ -V[3], -V[3], V[5], V[5] ] == I EE^2 / SW^2 *
{ 
 { 2 },
 { -1 },
 { -1 } 
},
  (*------    W-  Z  ~V+  ~V0  ------*)
   C[ -V[3], V[2], V[5], V[6] ] == - I CW EE^2 / SW^2 *
{ 
 { -1 },
 { 2 },
 { -1 } 
},
  (*------    Z  Z  ~V+  ~V-  ------*)
   C[ V[2], V[2], V[5], -V[5] ] == - I CW^2 EE^2 / SW^2 *
{ 
 { 2 },
 { -1 },
 { -1 } 
},
  (*------    ~V+  ~V+  ~V-  ~V-  ------*)
   C[ V[5], V[5], -V[5], -V[5] ] == I EE^2 / SW^2 *
{ 
 { 2 },
 { -1 },
 { -1 } 
},
  (*------    ~V+  ~V-  ~V0  ~V0  ------*)
   C[ V[5], -V[5], V[6], V[6] ] == - I EE^2 / SW^2 *
{ 
 { 2 },
 { -1 },
 { -1 } 
} 
}

M$LastModelRules = {}

Scan[ (RealQ[#] = True)&,
  { pi, alphaSMZ, EE, Q, Mm, Mtau, Ms, McMc, MbMb, Mtop, 
    MH, MZ, MW, wtop, wZ, wW, GG, Maux, Mcp, Mbp, 
    alphaE0, CW, SW, GF, vv, s12, s23, s13, c12, c23, 
    c13, LamQCD, Mb, Mt, Mc, a, dM, MV0, MVp, Mu, 
    Md, Me, Vud, Vus, Vub, Vcd, Vcs, Vcb, Vtd, Vts, 
    Vtb, lam  } ]





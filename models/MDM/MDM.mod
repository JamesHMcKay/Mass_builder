(* ----------------------------------------------------------------------------- *) 
(* This model file was automatically created by SARAH version4.7.0  *) 
(* SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223 *) 
(* (c) Florian Staub, 2013  *) 
(* ----------------------------------------------------------------------------- *) 
(* File created at 12:52 on 17.3.2017  *) 
(* ---------------------------------------------------------------------- *) 
 
 
IndexRange[  Index[Colour]  ] =NoUnfold[Range[3]]; 
IndexStyle[  Index[Colour, i_Integer ] ] := Greek[i];  
IndexRange[  Index[I3Gen]  ] =Range[3]; 
IndexStyle[  Index[I3Gen, i_Integer ] ] := Alph[ 8+i];  
IndexRange[  Index[Gluon]  ] =NoUnfold[Range[8]]; 
IndexStyle[  Index[Gluon, i_Integer ] ] := Alph[ 8+i];  

 
(* Definitions for trigonometric functions  
Sin[ThetaW]: STW
Sin[2*ThetaW]: S2TW
Cos[ThetaW]: CTW
Cos[2*ThetaW]: C2TW
*) 
 
Conjugate[STW] ^= STW
Conjugate[S2TW] ^= S2TW
Conjugate[CTW] ^= CTW
Conjugate[C2TW] ^= C2TW
 
 
Lam[a_,b_,c_]:=2*SUNT[a,b,c]; 
fSU3[a_,b_,c_]:=SUNF[a,b,c]; 
LambdaProd[a_,b_][c_,d_]:=4*SUNT[a,b,c,d]; 
 
 
M$ClassesDescription= {
S[3] == {SelfConjugate -> False,
Indices -> {},
Mass -> MassHp,
PropagatorLabel->ComposedChar["H","+"],
PropagatorType -> ScalarDash,
PropagatorArrow -> Forward},

 
S[2] == {SelfConjugate -> True,
Indices -> {},
Mass -> MassAh,
PropagatorLabel->ComposedChar["A","0"],
PropagatorType -> ScalarDash,
PropagatorArrow -> None},

 
S[1] == {SelfConjugate -> True,
Indices -> {},
Mass -> Masshh,
PropagatorLabel->h,
PropagatorType -> ScalarDash,
PropagatorArrow -> None},

 
F[1] == {SelfConjugate -> False,
Indices -> {Index[I3Gen]},
Mass -> MassFv,
PropagatorLabel->ComposedChar["\\nu",Index[I3Gen]],
PropagatorType -> Straight,
PropagatorArrow -> Forward},

 
F[5] == {SelfConjugate -> False,
Indices -> {},
Mass -> MChi,
PropagatorLabel->ComposedChar["\\chi","+"],
PropagatorType -> Straight,
PropagatorArrow -> Forward},

 
F[6] == {SelfConjugate -> False,
Indices -> {},
Mass -> MChi,
PropagatorLabel->ComposedChar["\\chi","++"],
PropagatorType -> Straight,
PropagatorArrow -> Forward},

 
F[7] == {SelfConjugate -> True,
Indices -> {},
Mass -> MChi,
PropagatorLabel->ComposedChar["\\chi","0"],
PropagatorType -> Straight,
PropagatorArrow -> None},

 
F[4] == {SelfConjugate -> False,
Indices -> {Index[I3Gen], Index[Colour]},
Mass -> MassFd,
PropagatorLabel->ComposedChar["d",Index[I3Gen],Index[Colour]],
PropagatorType -> Straight,
PropagatorArrow -> Forward},

 
F[3] == {SelfConjugate -> False,
Indices -> {Index[I3Gen], Index[Colour]},
Mass -> MassFu,
PropagatorLabel->ComposedChar["u",Index[I3Gen],Index[Colour]],
PropagatorType -> Straight,
PropagatorArrow -> Forward},

 
F[2] == {SelfConjugate -> False,
Indices -> {Index[I3Gen]},
Mass -> MassFe,
PropagatorLabel->ComposedChar["e",Index[I3Gen]],
PropagatorType -> Straight,
PropagatorArrow -> Forward},

 
V[5] == {SelfConjugate -> True,
Indices -> {Index[Gluon]},
Mass -> 0,
PropagatorLabel->ComposedChar["g",Index[Gluon]],
PropagatorType -> Sine,
PropagatorArrow -> None},

 
V[1] == {SelfConjugate -> True,
Indices -> {},
Mass -> 0,
PropagatorLabel->ComposedChar["\\gamma"],
PropagatorType -> Sine,
PropagatorArrow -> None},

 
V[2] == {SelfConjugate -> True,
Indices -> {},
Mass -> MassVZ,
PropagatorLabel->ComposedChar["Z"],
PropagatorType -> Sine,
PropagatorArrow -> None},

 
V[3] == {SelfConjugate -> False,
Indices -> {},
Mass -> MassVWp,
PropagatorLabel->ComposedChar["W","+"],
PropagatorType -> Sine,
PropagatorArrow -> Forward},

 
U[5] == {SelfConjugate -> False,
Indices -> {Index[Gluon]},
Mass -> 0,
PropagatorLabel->ComposedChar["\\eta",Index[Gluon],"G"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward},

 
U[1] == {SelfConjugate -> False,
Indices -> {},
Mass -> 0,
PropagatorLabel->ComposedChar["\\eta","\\gamma"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward},

 
U[2] == {SelfConjugate -> False,
Indices -> {},
Mass -> MassVZ,
PropagatorLabel->ComposedChar["\\eta","Z"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward},

 
U[3] == {SelfConjugate -> False,
Indices -> {},
Mass -> MassVWp,
PropagatorLabel->ComposedChar["\\eta","+"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward},

 
U[4] == {SelfConjugate -> False,
Indices -> {},
Mass -> MassVWp,
PropagatorLabel->ComposedChar["\\eta","-"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward}
 
}

 
MassFd[gen_, y_] = MassFd[gen]
MassFu[gen_, y_] = MassFu[gen]


GaugeXi[S[3,{a_Integer}]] = 1 /; a > 1 
GaugeXi[S[3,1]] = GaugeXi[VWp] /; a > 1 
GaugeXi[S[2,{a_Integer}]] = 1 /; a > 1 
GaugeXi[S[2,1]] = GaugeXi[VZ] /; a > 1 
GaugeXi[S[1,___]] = 1 


GaugeXi[V[5,___]] = GaugeXi[G]
GaugeXi[V[1,___]] = GaugeXi[P]
GaugeXi[V[2,___]] = GaugeXi[Z]
GaugeXi[V[3,___]] = GaugeXi[Wp]


M$CouplingMatrices= {
C[V[5, {ct1}], V[5, {ct2}], V[5, {ct3}], V[5, {ct4}]] == {{I*g3^2*(-(fSU3[1, ct1, ct4]*fSU3[1, ct2, ct3]) - fSU3[1, ct1, ct3]*fSU3[1, ct2, ct4] - fSU3[2, ct1, ct4]*fSU3[2, ct2, ct3] - fSU3[2, ct1, ct3]*fSU3[2, ct2, ct4] - fSU3[3, ct1, ct4]*fSU3[3, ct2, ct3] - fSU3[3, ct1, ct3]*fSU3[3, ct2, ct4] - fSU3[4, ct1, ct4]*fSU3[4, ct2, ct3] - fSU3[4, ct1, ct3]*fSU3[4, ct2, ct4] - fSU3[5, ct1, ct4]*fSU3[5, ct2, ct3] - fSU3[5, ct1, ct3]*fSU3[5, ct2, ct4] - fSU3[6, ct1, ct4]*fSU3[6, ct2, ct3] - fSU3[6, ct1, ct3]*fSU3[6, ct2, ct4] - fSU3[7, ct1, ct4]*fSU3[7, ct2, ct3] - fSU3[7, ct1, ct3]*fSU3[7, ct2, ct4] - fSU3[8, ct1, ct4]*fSU3[8, ct2, ct3] - fSU3[8, ct1, ct3]*fSU3[8, ct2, ct4])}, {I*g3^2*(fSU3[1, ct1, ct4]*fSU3[1, ct2, ct3] - fSU3[1, ct1, ct2]*fSU3[1, ct3, ct4] + fSU3[2, ct1, ct4]*fSU3[2, ct2, ct3] - fSU3[2, ct1, ct2]*fSU3[2, ct3, ct4] + fSU3[3, ct1, ct4]*fSU3[3, ct2, ct3] - fSU3[3, ct1, ct2]*fSU3[3, ct3, ct4] + fSU3[4, ct1, ct4]*fSU3[4, ct2, ct3] - fSU3[4, ct1, ct2]*fSU3[4, ct3, ct4] + fSU3[5, ct1, ct4]*fSU3[5, ct2, ct3] - fSU3[5, ct1, ct2]*fSU3[5, ct3, ct4] + fSU3[6, ct1, ct4]*fSU3[6, ct2, ct3] - fSU3[6, ct1, ct2]*fSU3[6, ct3, ct4] + fSU3[7, ct1, ct4]*fSU3[7, ct2, ct3] - fSU3[7, ct1, ct2]*fSU3[7, ct3, ct4] + fSU3[8, ct1, ct4]*fSU3[8, ct2, ct3] - fSU3[8, ct1, ct2]*fSU3[8, ct3, ct4])}, {I*g3^2*(fSU3[1, ct1, ct3]*fSU3[1, ct2, ct4] + fSU3[1, ct1, ct2]*fSU3[1, ct3, ct4] + fSU3[2, ct1, ct3]*fSU3[2, ct2, ct4] + fSU3[2, ct1, ct2]*fSU3[2, ct3, ct4] + fSU3[3, ct1, ct3]*fSU3[3, ct2, ct4] + fSU3[3, ct1, ct2]*fSU3[3, ct3, ct4] + fSU3[4, ct1, ct3]*fSU3[4, ct2, ct4] + fSU3[4, ct1, ct2]*fSU3[4, ct3, ct4] + fSU3[5, ct1, ct3]*fSU3[5, ct2, ct4] + fSU3[5, ct1, ct2]*fSU3[5, ct3, ct4] + fSU3[6, ct1, ct3]*fSU3[6, ct2, ct4] + fSU3[6, ct1, ct2]*fSU3[6, ct3, ct4] + fSU3[7, ct1, ct3]*fSU3[7, ct2, ct4] + fSU3[7, ct1, ct2]*fSU3[7, ct3, ct4] + fSU3[8, ct1, ct3]*fSU3[8, ct2, ct4] + fSU3[8, ct1, ct2]*fSU3[8, ct3, ct4])}},
 C[-V[3], V[1], V[1], V[3]] == {{I*g2^2*STW^2}, {I*g2^2*STW^2}, {(-2*I)*g2^2*STW^2}},
 C[-V[3], V[1], V[3], V[2]] == {{(I/2)*g2^2*S2TW}, {(-I)*g2^2*S2TW}, {(I/2)*g2^2*S2TW}},
 C[-V[3], -V[3], V[3], V[3]] == {{(2*I)*g2^2}, {(-I)*g2^2}, {(-I)*g2^2}},
 C[-V[3], V[3], V[2], V[2]] == {{(-2*I)*CTW^2*g2^2}, {I*CTW^2*g2^2}, {I*CTW^2*g2^2}},
 C[S[2], S[2], -V[3], V[3]] == {{(I/2)*g2^2}},
 C[S[2], S[2], V[2], V[2]] == {{(I/2)*CTW^2*g2^2 + I*CTW*g1*g2*STW + (I/2)*g1^2*STW^2}},
 C[S[2], S[3], -V[3], V[1]] == {{(CTW*g1*g2)/2}},
 C[S[2], S[3], -V[3], V[2]] == {{-(g1*g2*STW)/2}},
 C[S[2], -S[3], V[1], V[3]] == {{-(CTW*g1*g2)/2}},
 C[S[2], -S[3], V[3], V[2]] == {{(g1*g2*STW)/2}},
 C[S[1], S[1], -V[3], V[3]] == {{(I/2)*g2^2}},
 C[S[1], S[1], V[2], V[2]] == {{(I/2)*CTW^2*g2^2 + I*CTW*g1*g2*STW + (I/2)*g1^2*STW^2}},
 C[S[1], S[3], -V[3], V[1]] == {{(I/2)*CTW*g1*g2}},
 C[S[1], S[3], -V[3], V[2]] == {{(-I/2)*g1*g2*STW}},
 C[S[1], -S[3], V[1], V[3]] == {{(I/2)*CTW*g1*g2}},
 C[S[1], -S[3], V[3], V[2]] == {{(-I/2)*g1*g2*STW}},
 C[S[3], -S[3], V[1], V[1]] == {{(I/2)*CTW^2*g1^2 + I*CTW*g1*g2*STW + (I/2)*g2^2*STW^2}},
 C[S[3], -S[3], V[1], V[2]] == {{(I/2)*C2TW*g1*g2 - (I/4)*g1^2*S2TW + (I/4)*g2^2*S2TW}},
 C[S[3], -S[3], -V[3], V[3]] == {{(I/2)*g2^2}},
 C[S[3], -S[3], V[2], V[2]] == {{(I/2)*CTW^2*g2^2 - I*CTW*g1*g2*STW + (I/2)*g1^2*STW^2}},
 C[-F[4, {gt1, ct1}], F[4, {gt2, ct2}], S[2]] == {{-((IndexDelta[ct1, ct2]*(Conjugate[ZDL[gt2, 1]]*(Conjugate[ZDR[gt1, 1]]*Yd[1, 1] + Conjugate[ZDR[gt1, 2]]*Yd[2, 1] + Conjugate[ZDR[gt1, 3]]*Yd[3, 1]) + Conjugate[ZDL[gt2, 2]]*(Conjugate[ZDR[gt1, 1]]*Yd[1, 2] + Conjugate[ZDR[gt1, 2]]*Yd[2, 2] + Conjugate[ZDR[gt1, 3]]*Yd[3, 2]) + Conjugate[ZDL[gt2, 3]]*(Conjugate[ZDR[gt1, 1]]*Yd[1, 3] + Conjugate[ZDR[gt1, 2]]*Yd[2, 3] + Conjugate[ZDR[gt1, 3]]*Yd[3, 3])))/Sqrt[2])}, {(IndexDelta[ct1, ct2]*(ZDL[gt1, 1]*(Conjugate[Yd[1, 1]]*ZDR[gt2, 1] + Conjugate[Yd[2, 1]]*ZDR[gt2, 2] + Conjugate[Yd[3, 1]]*ZDR[gt2, 3]) + ZDL[gt1, 2]*(Conjugate[Yd[1, 2]]*ZDR[gt2, 1] + Conjugate[Yd[2, 2]]*ZDR[gt2, 2] + Conjugate[Yd[3, 2]]*ZDR[gt2, 3]) + ZDL[gt1, 3]*(Conjugate[Yd[1, 3]]*ZDR[gt2, 1] + Conjugate[Yd[2, 3]]*ZDR[gt2, 2] + Conjugate[Yd[3, 3]]*ZDR[gt2, 3])))/Sqrt[2]}},
 C[-F[2, {gt1}], F[2, {gt2}], S[2]] == {{-((Conjugate[ZEL[gt2, 1]]*(Conjugate[ZER[gt1, 1]]*Ye[1, 1] + Conjugate[ZER[gt1, 2]]*Ye[2, 1] + Conjugate[ZER[gt1, 3]]*Ye[3, 1]) + Conjugate[ZEL[gt2, 2]]*(Conjugate[ZER[gt1, 1]]*Ye[1, 2] + Conjugate[ZER[gt1, 2]]*Ye[2, 2] + Conjugate[ZER[gt1, 3]]*Ye[3, 2]) + Conjugate[ZEL[gt2, 3]]*(Conjugate[ZER[gt1, 1]]*Ye[1, 3] + Conjugate[ZER[gt1, 2]]*Ye[2, 3] + Conjugate[ZER[gt1, 3]]*Ye[3, 3]))/Sqrt[2])}, {(ZEL[gt1, 1]*(Conjugate[Ye[1, 1]]*ZER[gt2, 1] + Conjugate[Ye[2, 1]]*ZER[gt2, 2] + Conjugate[Ye[3, 1]]*ZER[gt2, 3]) + ZEL[gt1, 2]*(Conjugate[Ye[1, 2]]*ZER[gt2, 1] + Conjugate[Ye[2, 2]]*ZER[gt2, 2] + Conjugate[Ye[3, 2]]*ZER[gt2, 3]) + ZEL[gt1, 3]*(Conjugate[Ye[1, 3]]*ZER[gt2, 1] + Conjugate[Ye[2, 3]]*ZER[gt2, 2] + Conjugate[Ye[3, 3]]*ZER[gt2, 3]))/Sqrt[2]}},
 C[-F[3, {gt1, ct1}], F[3, {gt2, ct2}], S[2]] == {{-((IndexDelta[ct1, ct2]*(Conjugate[ZUL[gt2, 1]]*(Conjugate[ZUR[gt1, 1]]*Yu[1, 1] + Conjugate[ZUR[gt1, 2]]*Yu[2, 1] + Conjugate[ZUR[gt1, 3]]*Yu[3, 1]) + Conjugate[ZUL[gt2, 2]]*(Conjugate[ZUR[gt1, 1]]*Yu[1, 2] + Conjugate[ZUR[gt1, 2]]*Yu[2, 2] + Conjugate[ZUR[gt1, 3]]*Yu[3, 2]) + Conjugate[ZUL[gt2, 3]]*(Conjugate[ZUR[gt1, 1]]*Yu[1, 3] + Conjugate[ZUR[gt1, 2]]*Yu[2, 3] + Conjugate[ZUR[gt1, 3]]*Yu[3, 3])))/Sqrt[2])}, {(IndexDelta[ct1, ct2]*(ZUL[gt1, 1]*(Conjugate[Yu[1, 1]]*ZUR[gt2, 1] + Conjugate[Yu[2, 1]]*ZUR[gt2, 2] + Conjugate[Yu[3, 1]]*ZUR[gt2, 3]) + ZUL[gt1, 2]*(Conjugate[Yu[1, 2]]*ZUR[gt2, 1] + Conjugate[Yu[2, 2]]*ZUR[gt2, 2] + Conjugate[Yu[3, 2]]*ZUR[gt2, 3]) + ZUL[gt1, 3]*(Conjugate[Yu[1, 3]]*ZUR[gt2, 1] + Conjugate[Yu[2, 3]]*ZUR[gt2, 2] + Conjugate[Yu[3, 3]]*ZUR[gt2, 3])))/Sqrt[2]}},
 C[-F[4, {gt1, ct1}], F[4, {gt2, ct2}], S[1]] == {{((-I)*IndexDelta[ct1, ct2]*(Conjugate[ZDL[gt2, 1]]*(Conjugate[ZDR[gt1, 1]]*Yd[1, 1] + Conjugate[ZDR[gt1, 2]]*Yd[2, 1] + Conjugate[ZDR[gt1, 3]]*Yd[3, 1]) + Conjugate[ZDL[gt2, 2]]*(Conjugate[ZDR[gt1, 1]]*Yd[1, 2] + Conjugate[ZDR[gt1, 2]]*Yd[2, 2] + Conjugate[ZDR[gt1, 3]]*Yd[3, 2]) + Conjugate[ZDL[gt2, 3]]*(Conjugate[ZDR[gt1, 1]]*Yd[1, 3] + Conjugate[ZDR[gt1, 2]]*Yd[2, 3] + Conjugate[ZDR[gt1, 3]]*Yd[3, 3])))/Sqrt[2]}, {((-I)*IndexDelta[ct1, ct2]*(ZDL[gt1, 1]*(Conjugate[Yd[1, 1]]*ZDR[gt2, 1] + Conjugate[Yd[2, 1]]*ZDR[gt2, 2] + Conjugate[Yd[3, 1]]*ZDR[gt2, 3]) + ZDL[gt1, 2]*(Conjugate[Yd[1, 2]]*ZDR[gt2, 1] + Conjugate[Yd[2, 2]]*ZDR[gt2, 2] + Conjugate[Yd[3, 2]]*ZDR[gt2, 3]) + ZDL[gt1, 3]*(Conjugate[Yd[1, 3]]*ZDR[gt2, 1] + Conjugate[Yd[2, 3]]*ZDR[gt2, 2] + Conjugate[Yd[3, 3]]*ZDR[gt2, 3])))/Sqrt[2]}},
 C[-F[3, {gt1, ct1}], F[4, {gt2, ct2}], S[3]] == {{(-I)*IndexDelta[ct1, ct2]*(Conjugate[ZDL[gt2, 1]]*(Conjugate[ZUR[gt1, 1]]*Yu[1, 1] + Conjugate[ZUR[gt1, 2]]*Yu[2, 1] + Conjugate[ZUR[gt1, 3]]*Yu[3, 1]) + Conjugate[ZDL[gt2, 2]]*(Conjugate[ZUR[gt1, 1]]*Yu[1, 2] + Conjugate[ZUR[gt1, 2]]*Yu[2, 2] + Conjugate[ZUR[gt1, 3]]*Yu[3, 2]) + Conjugate[ZDL[gt2, 3]]*(Conjugate[ZUR[gt1, 1]]*Yu[1, 3] + Conjugate[ZUR[gt1, 2]]*Yu[2, 3] + Conjugate[ZUR[gt1, 3]]*Yu[3, 3]))}, {(-I)*IndexDelta[ct1, ct2]*((Conjugate[Yd[1, 1]]*ZDR[gt2, 1] + Conjugate[Yd[2, 1]]*ZDR[gt2, 2] + Conjugate[Yd[3, 1]]*ZDR[gt2, 3])*ZUL[gt1, 1] + (Conjugate[Yd[1, 2]]*ZDR[gt2, 1] + Conjugate[Yd[2, 2]]*ZDR[gt2, 2] + Conjugate[Yd[3, 2]]*ZDR[gt2, 3])*ZUL[gt1, 2] + (Conjugate[Yd[1, 3]]*ZDR[gt2, 1] + Conjugate[Yd[2, 3]]*ZDR[gt2, 2] + Conjugate[Yd[3, 3]]*ZDR[gt2, 3])*ZUL[gt1, 3])}},
 C[-F[2, {gt1}], F[2, {gt2}], S[1]] == {{((-I)*(Conjugate[ZEL[gt2, 1]]*(Conjugate[ZER[gt1, 1]]*Ye[1, 1] + Conjugate[ZER[gt1, 2]]*Ye[2, 1] + Conjugate[ZER[gt1, 3]]*Ye[3, 1]) + Conjugate[ZEL[gt2, 2]]*(Conjugate[ZER[gt1, 1]]*Ye[1, 2] + Conjugate[ZER[gt1, 2]]*Ye[2, 2] + Conjugate[ZER[gt1, 3]]*Ye[3, 2]) + Conjugate[ZEL[gt2, 3]]*(Conjugate[ZER[gt1, 1]]*Ye[1, 3] + Conjugate[ZER[gt1, 2]]*Ye[2, 3] + Conjugate[ZER[gt1, 3]]*Ye[3, 3])))/Sqrt[2]}, {((-I)*(ZEL[gt1, 1]*(Conjugate[Ye[1, 1]]*ZER[gt2, 1] + Conjugate[Ye[2, 1]]*ZER[gt2, 2] + Conjugate[Ye[3, 1]]*ZER[gt2, 3]) + ZEL[gt1, 2]*(Conjugate[Ye[1, 2]]*ZER[gt2, 1] + Conjugate[Ye[2, 2]]*ZER[gt2, 2] + Conjugate[Ye[3, 2]]*ZER[gt2, 3]) + ZEL[gt1, 3]*(Conjugate[Ye[1, 3]]*ZER[gt2, 1] + Conjugate[Ye[2, 3]]*ZER[gt2, 2] + Conjugate[Ye[3, 3]]*ZER[gt2, 3])))/Sqrt[2]}},
 C[-F[1, {gt1}], F[2, {gt2}], S[3]] == {{0}, {(-I)*(Conjugate[Ye[1, gt1]]*ZER[gt2, 1] + Conjugate[Ye[2, gt1]]*ZER[gt2, 2] + Conjugate[Ye[3, gt1]]*ZER[gt2, 3])}},
 C[-F[3, {gt1, ct1}], F[3, {gt2, ct2}], S[1]] == {{(I*IndexDelta[ct1, ct2]*(Conjugate[ZUL[gt2, 1]]*(Conjugate[ZUR[gt1, 1]]*Yu[1, 1] + Conjugate[ZUR[gt1, 2]]*Yu[2, 1] + Conjugate[ZUR[gt1, 3]]*Yu[3, 1]) + Conjugate[ZUL[gt2, 2]]*(Conjugate[ZUR[gt1, 1]]*Yu[1, 2] + Conjugate[ZUR[gt1, 2]]*Yu[2, 2] + Conjugate[ZUR[gt1, 3]]*Yu[3, 2]) + Conjugate[ZUL[gt2, 3]]*(Conjugate[ZUR[gt1, 1]]*Yu[1, 3] + Conjugate[ZUR[gt1, 2]]*Yu[2, 3] + Conjugate[ZUR[gt1, 3]]*Yu[3, 3])))/Sqrt[2]}, {(I*IndexDelta[ct1, ct2]*(ZUL[gt1, 1]*(Conjugate[Yu[1, 1]]*ZUR[gt2, 1] + Conjugate[Yu[2, 1]]*ZUR[gt2, 2] + Conjugate[Yu[3, 1]]*ZUR[gt2, 3]) + ZUL[gt1, 2]*(Conjugate[Yu[1, 2]]*ZUR[gt2, 1] + Conjugate[Yu[2, 2]]*ZUR[gt2, 2] + Conjugate[Yu[3, 2]]*ZUR[gt2, 3]) + ZUL[gt1, 3]*(Conjugate[Yu[1, 3]]*ZUR[gt2, 1] + Conjugate[Yu[2, 3]]*ZUR[gt2, 2] + Conjugate[Yu[3, 3]]*ZUR[gt2, 3])))/Sqrt[2]}},
 C[-F[4, {gt1, ct1}], F[3, {gt2, ct2}], -S[3]] == {{(-I)*IndexDelta[ct1, ct2]*(Conjugate[ZUL[gt2, 1]]*(Conjugate[ZDR[gt1, 1]]*Yd[1, 1] + Conjugate[ZDR[gt1, 2]]*Yd[2, 1] + Conjugate[ZDR[gt1, 3]]*Yd[3, 1]) + Conjugate[ZUL[gt2, 2]]*(Conjugate[ZDR[gt1, 1]]*Yd[1, 2] + Conjugate[ZDR[gt1, 2]]*Yd[2, 2] + Conjugate[ZDR[gt1, 3]]*Yd[3, 2]) + Conjugate[ZUL[gt2, 3]]*(Conjugate[ZDR[gt1, 1]]*Yd[1, 3] + Conjugate[ZDR[gt1, 2]]*Yd[2, 3] + Conjugate[ZDR[gt1, 3]]*Yd[3, 3]))}, {(-I)*IndexDelta[ct1, ct2]*(ZDL[gt1, 1]*(Conjugate[Yu[1, 1]]*ZUR[gt2, 1] + Conjugate[Yu[2, 1]]*ZUR[gt2, 2] + Conjugate[Yu[3, 1]]*ZUR[gt2, 3]) + ZDL[gt1, 2]*(Conjugate[Yu[1, 2]]*ZUR[gt2, 1] + Conjugate[Yu[2, 2]]*ZUR[gt2, 2] + Conjugate[Yu[3, 2]]*ZUR[gt2, 3]) + ZDL[gt1, 3]*(Conjugate[Yu[1, 3]]*ZUR[gt2, 1] + Conjugate[Yu[2, 3]]*ZUR[gt2, 2] + Conjugate[Yu[3, 3]]*ZUR[gt2, 3]))}},
 C[-F[2, {gt1}], F[1, {gt2}], -S[3]] == {{(-I)*(Conjugate[ZER[gt1, 1]]*Ye[1, gt2] + Conjugate[ZER[gt1, 2]]*Ye[2, gt2] + Conjugate[ZER[gt1, 3]]*Ye[3, gt2])}, {0}},
 C[F[7], F[5], -V[3]] == {{I*Sqrt[3]*g2}, {I*Sqrt[3]*g2}},
 C[-F[5], F[5], V[1]] == {{(-I)*g2*STW}, {(-I)*g2*STW}},
 C[-F[6], F[5], V[3]] == {{I*Sqrt[2]*g2}, {I*Sqrt[2]*g2}},
 C[-F[5], F[5], V[2]] == {{(-I)*CTW*g2}, {(-I)*CTW*g2}},
 C[-F[4, {gt1, ct1}], F[4, {gt2, ct2}], V[5, {ct3}]] == {{(-I/2)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct1, ct2]}, {(-I/2)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct1, ct2]}},
 C[-F[4, {gt1, ct1}], F[4, {gt2, ct2}], V[1]] == {{(-I/6)*(CTW*g1 - 3*g2*STW)*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}, {(I/3)*CTW*g1*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}},
 C[-F[3, {gt1, ct1}], F[4, {gt2, ct2}], V[3]] == {{((-I)*g2*IndexDelta[ct1, ct2]*(Conjugate[ZDL[gt2, 1]]*ZUL[gt1, 1] + Conjugate[ZDL[gt2, 2]]*ZUL[gt1, 2] + Conjugate[ZDL[gt2, 3]]*ZUL[gt1, 3]))/Sqrt[2]}, {0}},
 C[-F[4, {gt1, ct1}], F[4, {gt2, ct2}], V[2]] == {{(I/6)*(3*CTW*g2 + g1*STW)*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}, {(-I/3)*g1*STW*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}},
 C[-F[2, {gt1}], F[2, {gt2}], V[1]] == {{(I/2)*(CTW*g1 + g2*STW)*IndexDelta[gt1, gt2]}, {I*CTW*g1*IndexDelta[gt1, gt2]}},
 C[-F[1, {gt1}], F[2, {gt2}], V[3]] == {{((-I)*g2*Conjugate[ZEL[gt2, gt1]])/Sqrt[2]}, {0}},
 C[-F[2, {gt1}], F[2, {gt2}], V[2]] == {{(I/2)*(CTW*g2 - g1*STW)*IndexDelta[gt1, gt2]}, {(-I)*g1*STW*IndexDelta[gt1, gt2]}},
 C[-F[6], F[6], V[1]] == {{(-2*I)*g2*STW}, {(-2*I)*g2*STW}},
 C[-F[6], F[6], V[2]] == {{(-2*I)*CTW*g2}, {(-2*I)*CTW*g2}},
 C[-F[5], F[6], -V[3]] == {{I*Sqrt[2]*g2}, {I*Sqrt[2]*g2}},
 C[-F[5], F[7], V[3]] == {{I*Sqrt[3]*g2}, {I*Sqrt[3]*g2}},
 C[-F[3, {gt1, ct1}], F[3, {gt2, ct2}], V[5, {ct3}]] == {{(-I/2)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct1, ct2]}, {(-I/2)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct1, ct2]}},
 C[-F[3, {gt1, ct1}], F[3, {gt2, ct2}], V[1]] == {{(-I/6)*(CTW*g1 + 3*g2*STW)*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}, {((-2*I)/3)*CTW*g1*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}},
 C[-F[3, {gt1, ct1}], F[3, {gt2, ct2}], V[2]] == {{(-I/6)*(3*CTW*g2 - g1*STW)*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}, {((2*I)/3)*g1*STW*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}},
 C[-F[4, {gt1, ct1}], F[3, {gt2, ct2}], -V[3]] == {{((-I)*g2*IndexDelta[ct1, ct2]*(Conjugate[ZUL[gt2, 1]]*ZDL[gt1, 1] + Conjugate[ZUL[gt2, 2]]*ZDL[gt1, 2] + Conjugate[ZUL[gt2, 3]]*ZDL[gt1, 3]))/Sqrt[2]}, {0}},
 C[-F[1, {gt1}], F[1, {gt2}], V[2]] == {{(-I/2)*(CTW*g2 + g1*STW)*IndexDelta[gt1, gt2]}, {0}},
 C[-F[2, {gt1}], F[1, {gt2}], -V[3]] == {{((-I)*g2*ZEL[gt1, gt2])/Sqrt[2]}, {0}},
 C[S[2], S[2], S[1]] == {{(-2*I)*Lam*v}},
 C[S[1], S[1], S[1]] == {{(-6*I)*Lam*v}},
 C[S[1], S[3], -S[3]] == {{(-2*I)*Lam*v}},
 C[S[2], S[1], V[2]] == {{(-(CTW*g2) - g1*STW)/2}},
 C[S[2], S[3], -V[3]] == {{g2/2}},
 C[S[2], -S[3], V[3]] == {{g2/2}},
 C[S[1], S[3], -V[3]] == {{(I/2)*g2}},
 C[S[1], -S[3], V[3]] == {{(-I/2)*g2}},
 C[S[3], -S[3], V[1]] == {{(-I/2)*(CTW*g1 + g2*STW)}},
 C[S[3], -S[3], V[2]] == {{(-I/2)*(CTW*g2 - g1*STW)}},
 C[S[1], -V[3], V[3]] == {{(I/2)*g2^2*v}},
 C[S[1], V[2], V[2]] == {{(I/2)*(CTW*g2 + g1*STW)^2*v}},
 C[S[3], -V[3], V[1]] == {{(I/2)*CTW*g1*g2*v}},
 C[S[3], -V[3], V[2]] == {{(-I/2)*g1*g2*STW*v}},
 C[-S[3], V[1], V[3]] == {{(I/2)*CTW*g1*g2*v}},
 C[-S[3], V[3], V[2]] == {{(-I/2)*g1*g2*STW*v}},
 C[S[2], S[2], S[2], S[2]] == {{(-6*I)*Lam}},
 C[S[2], S[2], S[1], S[1]] == {{(-2*I)*Lam}},
 C[S[2], S[2], S[3], -S[3]] == {{(-2*I)*Lam}},
 C[S[1], S[1], S[1], S[1]] == {{(-6*I)*Lam}},
 C[S[1], S[1], S[3], -S[3]] == {{(-2*I)*Lam}},
 C[S[3], S[3], -S[3], -S[3]] == {{(-4*I)*Lam}},
 C[V[5, {ct1}], V[5, {ct2}], V[5, {ct3}]] == {{g3*fSU3[ct1, ct2, ct3]}},
 C[-V[3], V[1], V[3]] == {{(-I)*g2*STW}},
 C[-V[3], V[3], V[2]] == {{I*CTW*g2}},
 C[S[2], -U[3], U[3]] == {{(g2^2*v*GaugeXi[Wp])/4}},
 C[S[2], -U[4], U[4]] == {{-(g2^2*v*GaugeXi[Wp])/4}},
 C[S[1], -U[2], U[1]] == {{(I/8)*(2*C2TW*g1*g2 + (g1^2 - g2^2)*S2TW)*v*GaugeXi[Z]}},
 C[S[3], -U[3], U[1]] == {{(-I/4)*g2*(CTW*g1 + g2*STW)*v*GaugeXi[Wp]}},
 C[-S[3], -U[4], U[1]] == {{(-I/4)*g2*(CTW*g1 + g2*STW)*v*GaugeXi[Wp]}},
 C[S[1], -U[3], U[3]] == {{(-I/4)*g2^2*v*GaugeXi[Wp]}},
 C[-S[3], -U[2], U[3]] == {{(I/4)*g2*(CTW*g2 + g1*STW)*v*GaugeXi[Z]}},
 C[S[1], -U[4], U[4]] == {{(-I/4)*g2^2*v*GaugeXi[Wp]}},
 C[S[3], -U[2], U[4]] == {{(I/4)*g2*(CTW*g2 + g1*STW)*v*GaugeXi[Z]}},
 C[S[1], -U[2], U[2]] == {{(-I/4)*(CTW*g2 + g1*STW)^2*v*GaugeXi[Z]}},
 C[S[3], -U[3], U[2]] == {{(-I/4)*g2*(CTW*g2 - g1*STW)*v*GaugeXi[Wp]}},
 C[-S[3], -U[4], U[2]] == {{(-I/4)*g2*(CTW*g2 - g1*STW)*v*GaugeXi[Wp]}},
 C[-U[5, {ct1}], U[5, {ct2}], V[5, {ct3}]] == {{g3*fSU3[ct1, ct2, ct3]}, {0}},
 C[-U[3], U[1], V[3]] == {{(-I)*g2*STW}, {0}},
 C[-U[4], U[1], -V[3]] == {{I*g2*STW}, {0}},
 C[-U[3], U[3], V[1]] == {{I*g2*STW}, {0}},
 C[-U[3], U[3], V[2]] == {{I*CTW*g2}, {0}},
 C[-U[1], U[3], -V[3]] == {{(-I)*g2*STW}, {0}},
 C[-U[2], U[3], -V[3]] == {{(-I)*CTW*g2}, {0}},
 C[-U[4], U[4], V[1]] == {{(-I)*g2*STW}, {0}},
 C[-U[1], U[4], V[3]] == {{I*g2*STW}, {0}},
 C[-U[2], U[4], V[3]] == {{I*CTW*g2}, {0}},
 C[-U[4], U[4], V[2]] == {{(-I)*CTW*g2}, {0}},
 C[-U[3], U[2], V[3]] == {{(-I)*CTW*g2}, {0}},
 C[-U[4], U[2], -V[3]] == {{I*CTW*g2}, {0}}
 }

 
Conjugate[g1] ^= g1; 
Conjugate[g2] ^= g2; 
Conjugate[g3] ^= g3; 
Conjugate[v] ^= v; 
Conjugate[ZZ[a___]] ^= ZZ[a]; 
Conjugate[aEWinv] ^= aEWinv; 
Conjugate[Gf] ^= Gf; 

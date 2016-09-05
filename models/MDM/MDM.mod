(* ----------------------------------------------------------------------------- *) 
(* This model file was automatically created by SARAH version4.7.0  *) 
(* SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223 *) 
(* (c) Florian Staub, 2013  *) 
(* ----------------------------------------------------------------------------- *) 
(* File created at 13:03 on 9.5.2016  *) 
(* ---------------------------------------------------------------------- *) 
 
 
IndexRange[  Index[Colour]  ] =NoUnfold[Range[3]];
IndexStyle[  Index[Colour, i_Integer ] ] := Greek[i];  
IndexRange[  Index[I3Gen]  ] =Range[3]; 
IndexStyle[  Index[I3Gen, i_Integer ] ] := Alph[ 8+i];  
IndexRange[  Index[Gluon]  ] =NoUnfold[Range[8]]; 
IndexStyle[  Index[Gluon, i_Integer ] ] := Alph[ 8+i];  

 
(* Definitions for trigonometric functions  
Sin[ThetaW]: sw
Sin[2*ThetaW]: sw2
Cos[ThetaW]: cw
Cos[2*ThetaW]: cw2
*) 
 
Conjugate[sw] ^= sw
Conjugate[sw2] ^= sw2
Conjugate[cw] ^= cw
Conjugate[cw2] ^= cw2
 
 
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
Mass -> mh,
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

 
F[6] == {SelfConjugate -> True,
Indices -> {},
Mass -> MChi,
PropagatorLabel->ComposedChar["\\chi","0"],
PropagatorType -> Straight,
PropagatorArrow -> None},

 
V[5] == {SelfConjugate -> True,
Indices -> {Index[Gluon]},
Mass -> mg,
PropagatorLabel->ComposedChar["g",Index[Gluon]],
PropagatorType -> Sine,
PropagatorArrow -> None},

 
V[1] == {SelfConjugate -> True,
Indices -> {},
Mass -> ma,
PropagatorLabel->ComposedChar["\\gamma"],
PropagatorType -> Sine,
PropagatorArrow -> None},

 
V[2] == {SelfConjugate -> True,
Indices -> {},
Mass -> mz,
PropagatorLabel->ComposedChar["Z"],
PropagatorType -> Sine,
PropagatorArrow -> None},

 
V[3] == {SelfConjugate -> False,
Indices -> {},
Mass -> mw,
PropagatorLabel->ComposedChar["W","+"],
PropagatorType -> Sine,
PropagatorArrow -> Forward},

 
U[5] == {SelfConjugate -> False,
Indices -> {Index[Gluon]},
Mass -> mg,
PropagatorLabel->ComposedChar["\\eta",Index[Gluon],"G"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward},

 
U[1] == {SelfConjugate -> False,
Indices -> {},
Mass -> ma,
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
 C[-V[3], V[1], V[1], V[3]] == {{I*g2^2*sw2}, {I*g2^2*sw2}, {(-2*I)*g2^2*sw2}},
 C[-V[3], V[1], V[3], V[2]] == {{(I/2)*g2^2*S2TW}, {(-I)*g2^2*S2TW}, {(I/2)*g2^2*S2TW}},
 C[-V[3], -V[3], V[3], V[3]] == {{(2*I)*g2^2}, {(-I)*g2^2}, {(-I)*g2^2}},
 C[-V[3], V[3], V[2], V[2]] == {{(-2*I)*cw2*g2^2}, {I*cw2*g2^2}, {I*cw2*g2^2}},
 C[S[2], S[2], -V[3], V[3]] == {{(I/2)*g2^2}},
 C[S[2], S[2], V[2], V[2]] == {{(I/2)*cw2*g2^2 + I*cw*g1*g2*sw + (I/2)*g1^2*sw2}},
 C[S[2], S[3], -V[3], V[1]] == {{(cw*g1*g2)/2}},
 C[S[2], S[3], -V[3], V[2]] == {{-(g1*g2*sw)/2}},
 C[S[2], -S[3], V[1], V[3]] == {{-(cw*g1*g2)/2}},
 C[S[2], -S[3], V[3], V[2]] == {{(g1*g2*sw)/2}},
 C[S[1], S[1], -V[3], V[3]] == {{(I/2)*g2^2}},
 C[S[1], S[1], V[2], V[2]] == {{(I/2)*cw2*g2^2 + I*cw*g1*g2*sw + (I/2)*g1^2*sw2}},
 C[S[1], S[3], -V[3], V[1]] == {{(I/2)*cw*g1*g2}},
 C[S[1], S[3], -V[3], V[2]] == {{(-I/2)*g1*g2*sw}},
 C[S[1], -S[3], V[1], V[3]] == {{(I/2)*cw*g1*g2}},
 C[S[1], -S[3], V[3], V[2]] == {{(-I/2)*g1*g2*sw}},
 C[S[3], -S[3], V[1], V[1]] == {{(I/2)*cw2*g1^2 + I*cw*g1*g2*sw + (I/2)*g2^2*sw2}},
 C[S[3], -S[3], V[1], V[2]] == {{(I/2)*C2TW*g1*g2 - (I/4)*g1^2*S2TW + (I/4)*g2^2*S2TW}},
 C[S[3], -S[3], -V[3], V[3]] == {{(I/2)*g2^2}},
 C[S[3], -S[3], V[2], V[2]] == {{(I/2)*cw2*g2^2 - I*cw*g1*g2*sw + (I/2)*g1^2*sw2}},
 C[F[6], F[5], -V[3]] == {{I*g2}, {I*g2}},
 C[-F[5], F[5], V[1]] == {{(-I)*g2*sw}, {(-I)*g2*sw}},
 C[-F[5], F[5], V[2]] == {{(-I)*cw*g2}, {(-I)*cw*g2}},
 C[-F[5], F[6], V[3]] == {{I*g2}, {I*g2}},
 C[S[2], S[2], S[1]] == {{(-2*I)*Lam*v}},
 C[S[1], S[1], S[1]] == {{(-6*I)*Lam*v}},
 C[S[1], S[3], -S[3]] == {{(-2*I)*Lam*v}},
 C[S[2], S[1], V[2]] == {{(-(cw*g2) - g1*sw)/2}},
 C[S[2], S[3], -V[3]] == {{g2/2}},
 C[S[2], -S[3], V[3]] == {{g2/2}},
 C[S[1], S[3], -V[3]] == {{(I/2)*g2}},
 C[S[1], -S[3], V[3]] == {{(-I/2)*g2}},
 C[S[3], -S[3], V[1]] == {{(-I/2)*(cw*g1 + g2*sw)}},
 C[S[3], -S[3], V[2]] == {{(-I/2)*(cw*g2 - g1*sw)}},
 C[S[1], -V[3], V[3]] == {{(I/2)*g2^2*v}},
 C[S[1], V[2], V[2]] == {{(I/2)*(cw*g2 + g1*sw)^2*v}},
 C[S[3], -V[3], V[1]] == {{(I/2)*cw*g1*g2*v}},
 C[S[3], -V[3], V[2]] == {{(-I/2)*g1*g2*sw*v}},
 C[-S[3], V[1], V[3]] == {{(I/2)*cw*g1*g2*v}},
 C[-S[3], V[3], V[2]] == {{(-I/2)*g1*g2*sw*v}},
 C[S[2], S[2], S[2], S[2]] == {{(-6*I)*Lam}},
 C[S[2], S[2], S[1], S[1]] == {{(-2*I)*Lam}},
 C[S[2], S[2], S[3], -S[3]] == {{(-2*I)*Lam}},
 C[S[1], S[1], S[1], S[1]] == {{(-6*I)*Lam}},
 C[S[1], S[1], S[3], -S[3]] == {{(-2*I)*Lam}},
 C[S[3], S[3], -S[3], -S[3]] == {{(-4*I)*Lam}},
 C[V[5, {ct1}], V[5, {ct2}], V[5, {ct3}]] == {{g3*fSU3[ct1, ct2, ct3]}},
 C[-V[3], V[1], V[3]] == {{(-I)*g2*sw}},
 C[-V[3], V[3], V[2]] == {{I*cw*g2}},
 C[S[2], -U[3], U[3]] == {{(g2^2*v*GaugeXi[Wp])/4}},
 C[S[2], -U[4], U[4]] == {{-(g2^2*v*GaugeXi[Wp])/4}},
 C[S[1], -U[2], U[1]] == {{(I/8)*(2*C2TW*g1*g2 + (g1^2 - g2^2)*S2TW)*v*GaugeXi[Z]}},
 C[S[3], -U[3], U[1]] == {{(-I/4)*g2*(cw*g1 + g2*sw)*v*GaugeXi[Wp]}},
 C[-S[3], -U[4], U[1]] == {{(-I/4)*g2*(cw*g1 + g2*sw)*v*GaugeXi[Wp]}},
 C[S[1], -U[3], U[3]] == {{(-I/4)*g2^2*v*GaugeXi[Wp]}},
 C[-S[3], -U[2], U[3]] == {{(I/4)*g2*(cw*g2 + g1*sw)*v*GaugeXi[Z]}},
 C[S[1], -U[4], U[4]] == {{(-I/4)*g2^2*v*GaugeXi[Wp]}},
 C[S[3], -U[2], U[4]] == {{(I/4)*g2*(cw*g2 + g1*sw)*v*GaugeXi[Z]}},
 C[S[1], -U[2], U[2]] == {{(-I/4)*(cw*g2 + g1*sw)^2*v*GaugeXi[Z]}},
 C[S[3], -U[3], U[2]] == {{(-I/4)*g2*(cw*g2 - g1*sw)*v*GaugeXi[Wp]}},
 C[-S[3], -U[4], U[2]] == {{(-I/4)*g2*(cw*g2 - g1*sw)*v*GaugeXi[Wp]}},
 C[-U[5, {ct1}], U[5, {ct2}], V[5, {ct3}]] == {{g3*fSU3[ct1, ct2, ct3]}, {0}},
 C[-U[3], U[1], V[3]] == {{(-I)*g2*sw}, {0}},
 C[-U[4], U[1], -V[3]] == {{I*g2*sw}, {0}},
 C[-U[3], U[3], V[1]] == {{I*g2*sw}, {0}},
 C[-U[3], U[3], V[2]] == {{I*cw*g2}, {0}},
 C[-U[1], U[3], -V[3]] == {{(-I)*g2*sw}, {0}},
 C[-U[2], U[3], -V[3]] == {{(-I)*cw*g2}, {0}},
 C[-U[4], U[4], V[1]] == {{(-I)*g2*sw}, {0}},
 C[-U[1], U[4], V[3]] == {{I*g2*sw}, {0}},
 C[-U[2], U[4], V[3]] == {{I*cw*g2}, {0}},
 C[-U[4], U[4], V[2]] == {{(-I)*cw*g2}, {0}},
 C[-U[3], U[2], V[3]] == {{(-I)*cw*g2}, {0}},
 C[-U[4], U[2], -V[3]] == {{I*cw*g2}, {0}}
 }

 
Conjugate[g1] ^= g1; 
Conjugate[g2] ^= g2; 
Conjugate[g3] ^= g3; 
Conjugate[v] ^= v; 
Conjugate[ZZ[a___]] ^= ZZ[a]; 
Conjugate[aEWinv] ^= aEWinv; 
Conjugate[Gf] ^= Gf; 

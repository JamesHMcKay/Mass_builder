(* ----------------------------------------------------------------------------- *)
(* This model file was automatically created by SARAH version4.7.0  *) 
(* SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223 *) 
(* (c) Florian Staub, 2013  *) 
(* ----------------------------------------------------------------------------- *) 
(* File created at 15:06 on 9.10.2016  *) 
(* ---------------------------------------------------------------------- *) 
 
 

 
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
Mass -> mw,
PropagatorLabel->ComposedChar["H","+"],
PropagatorType -> ScalarDash,
PropagatorArrow -> Forward},

 
S[2] == {SelfConjugate -> True,
Indices -> {},
Mass -> mz,
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
Indices -> {},
Mass -> MChi,
PropagatorLabel->ComposedChar["\\chi","+"],
PropagatorType -> Straight,
PropagatorArrow -> Forward},

 
F[2] == {SelfConjugate -> True,
Indices -> {},
Mass -> MChi,
PropagatorLabel->ComposedChar["\\chi","0"],
PropagatorType -> Straight,
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

 
U[1] == {SelfConjugate -> False,
Indices -> {},
Mass -> ma,
PropagatorLabel->ComposedChar["\\eta","\\gamma"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward},

 
U[2] == {SelfConjugate -> False,
Indices -> {},
Mass -> mz,
PropagatorLabel->ComposedChar["\\eta","Z"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward},

 
U[3] == {SelfConjugate -> False,
Indices -> {},
Mass -> mw,
PropagatorLabel->ComposedChar["\\eta","+"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward},

 
U[4] == {SelfConjugate -> False,
Indices -> {},
Mass -> mw,
PropagatorLabel->ComposedChar["\\eta","-"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward}
 
}

 


GaugeXi[S[3,{a_Integer}]] = 1 /; a > 1 
GaugeXi[S[3,1]] = GaugeXi[VWp] /; a > 1 
GaugeXi[S[2,{a_Integer}]] = 1 /; a > 1 
GaugeXi[S[2,1]] = GaugeXi[VZ] /; a > 1 
GaugeXi[S[1,___]] = 1 


GaugeXi[V[1,___]] = GaugeXi[P]
GaugeXi[V[2,___]] = GaugeXi[Z]
GaugeXi[V[3,___]] = GaugeXi[Wp]


M$CouplingMatrices= {
C[ -V[3], V[3] ] == I *
    { {0, dwZ, 0},
      {0, dwM, 0},
      {0, 0, 0} },
  C[ V[2], V[2] ] == I *
    { {0, dzZ, 0},
      {0, dzM, 0},
      {0, 0, 0} },
  C[ V[1], V[1] ] == I *
    { {0, dgammaZ, 0},
      {0, 0, 0},
      {0, 0,0} },
C[ V[1], V[2] ] == I *
    { {0, dZgammaZ,0},
      {0, dZgammaM,0},
      {0, 0,0}},
C[ F[1], -F[1]] == -I *  { {0, 0},{0, d1Z},{0,d1m},{0,0} },
C[ F[2], F[2]] == -I *  { {0, 0},{0, d2Z},{0,d2m},{0,0} },
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
 C[F[2], F[1], -V[3]] == {{I*g2,0}, {I*g2,I*dg2}},
 C[-F[1], F[1], V[1]] == {{(-I)*g2*STW}, {(-I)*g2*STW}},
 C[-F[1], F[1], V[2]] == {{(-I)*CTW*g2}, {(-I)*CTW*g2}},
 C[-F[1], F[2], V[3]] == {{I*g2,0}, {I*g2,I*dg2}},
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
Conjugate[v] ^= v; 
Conjugate[ZZ[a___]] ^= ZZ[a]; 
Conjugate[aEWinv] ^= aEWinv; 
Conjugate[Gf] ^= Gf; 

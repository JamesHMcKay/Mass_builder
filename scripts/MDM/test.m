(* ::Package:: *)

Quit[]


$LoadTARCER = True;
$LoadFeynArts = True;
<< FeynCalc/FeynCalc.m;
AppendTo[$Path, 
  "/Users/jamesmckay/Documents/Programs/Mass_builder/src/"];
<< MassBuilder.m;
SetOptions[DiracSlash, Dimension -> D, 
 FeynCalcInternal -> True]; SetOptions[DiracTrace, 
 DiracTraceEvaluate -> True]; null = 0;


(* Determine one-loop counter-term couplings *)

path = "/Users/jamesmckay/Documents/Programs/Mass_builder/";
kappa = 1/(16 \[Pi]^2);
sw=Sin[\[Theta]];
cw=Cos[\[Theta]];
STW=Sin[\[Theta]];
CTW=Cos[\[Theta]];

ClearScalarProducts[];

Get[FileNameJoin[{path, 
    "/models/MDM///output/math_data_F7_{1,2}_1.mx"}]];
SEMBtemp =  SelfEnergyFinite*kappa;


SEMB=FullSimplify[makeFiniteAmplitude[SEMBtemp,0,D]/.Pair[Momentum[p], Momentum[p]]->p^2]


SEFS = FullSimplify[kappa*(-8*MChi*(-(1/2) + B0[p^2,MChi^2,mw^2])*3*g2^2-4*p*3*g2^2*((1/2) + B1[p^2,MChi^2,mw^2]))]


SEFS2=FullSimplify[SEFS/.B0[p_,m1_,m2_]->  b*I*TBI[4, p, {{1, Sqrt[m1]}, {1, Sqrt[m2]}}]/.A0[m_]->a*I*TAI[4, 0, {{1, Sqrt[m]}}],{mw>0,MChi>0}]


FullSimplify[SEFS2-SEMB/.p->MChi/.a->-1/.b->-1]

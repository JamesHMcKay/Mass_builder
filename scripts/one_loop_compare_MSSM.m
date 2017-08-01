(* ::Package:: *)

Quit[]


$LoadTARCER = True;
$LoadFeynArts = True;
<< FeynCalc/FeynCalc.m;
AppendTo[$Path,"/Users/jamesmckay/Documents/Programs/Mass_builder/src/"];
<< MassBuilder.m;
SetOptions[DiracSlash, Dimension -> D, 
 FeynCalcInternal -> True]; SetOptions[DiracTrace, DiracTraceEvaluate -> True]; null = 0;


(* Determine one-loop counter-term couplings *)

path = "/Users/jamesmckay/Documents/Programs/Mass_builder_results/";
kappa = 1/(16 \[Pi]^2);
sw=Sin[\[Theta]];
cw=Cos[\[Theta]];
(*F11*)

F11SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_1_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_2_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_3_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

F11SEdiv = makeFiniteAmplitude[F11SE, -1, D];

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_1_1c.mx"}]];

F11ct = makeFiniteAmplitude[SelfEnergyFinite*kappa, -1, D] /. 
    MassBuilderCTM1 -> 0 /. MassBuilderCTZ1 -> 0;

eq1 = FullSimplify[Coefficient[F11ct + F11SEdiv, p]];
eq2 = FullSimplify[Coefficient[F11ct + F11SEdiv, p, 0]];
solF11 = Solve[{eq1 == 0, eq2 == 0}, {d2Z, d2m}]

(*F12*)

F12SE = 0;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F12_g1_1_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F12_g1_2_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F12_g1_3_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

F12SEdiv = makeFiniteAmplitude[F12SE, -1, D];
Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F12_g1_1_1c.mx"}]];
F12ct = makeFiniteAmplitude[SelfEnergyFinite*kappa, -1, D] /. 
    MassBuilderCTM1 -> 0 /. MassBuilderCTZ1 -> 0;
eq1 = FullSimplify[Coefficient[F12ct + F12SEdiv, p]];
eq2 = FullSimplify[Coefficient[F12ct + F12SEdiv, p, 0]];
solF12 = Solve[{eq1 == 0, eq2 == 0}, {d1Z, d1m}]
Set @@@ solF11[[1]];
Set @@@ solF12[[1]];


F11SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_1_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_2_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_3_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/massiveMSSM/output/math_data_F11_g1_1_1c.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;
F11SE = F11SE /. MassBuilderCTZ1 -> 0 /. MassBuilderCTZ2 -> 0 /. 
     MassBuilderCTM1 -> 0 /. MassBuilderCTM2 -> 0 /. 
   D -> 4 - 2 MassBuilderEpsilon;

C0 = Coefficient[FullSimplify[F11SE], MassBuilderEpsilon, 0];
C1 = Coefficient[FullSimplify[F11SE], MassBuilderEpsilon, 1];
C2 = Coefficient[FullSimplify[F11SE], MassBuilderEpsilon, 2];
SEtest = C0 + C1*MassBuilderEpsilon + C2*MassBuilderEpsilon^2;
FullSimplify[SEtest - F11SE];
SEN = FullSimplify[SEtest] /. MassBuilderEpsilon -> 0 /. Pair[Momentum[p], Momentum[p]] -> p^2;

(* express in terms of B1 *)
SEN2 = SEN /. TBI[4, p^2, {{1, mw}, {1, MChi}}] -> b Bwc[s]\
/. TAI[4, 0, {{1, mw}}] -> Aw /. TAI[4, 0, {{1, MChi}}] -> Ac\
/. TAI[4, 0, {{1, ma}}] -> Aa /. TAI[4, 0, {{1, mz}}] -> Az;
SEN2 = FullSimplify[SEN2 /. (Aw) -> 2 p^2 a B1cw[s] -b(p^2 + MChi^2 - mw^2)*Bwc[s] + Ac]



(* F12 *)
F12SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_1_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_2_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_3_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_1_1c.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

F12SE = F12SE/.MassBuilderCTZ1->0/.MassBuilderCTZ2->0/.MassBuilderCTM1->0/.MassBuilderCTM2->0/.D->4-2MassBuilderEpsilon;
C0=Coefficient[Simplify[F12SE],MassBuilderEpsilon,0];
C1=Coefficient[Simplify[F12SE],MassBuilderEpsilon,1];
C2=Coefficient[Simplify[F12SE],MassBuilderEpsilon,2];
SEtest2 = C0+C1*MassBuilderEpsilon+C2*MassBuilderEpsilon^2;
Simplify[SEtest2-F12SE];
SEC=Simplify[SEtest2]/.MassBuilderEpsilon->0/.Pair[Momentum[p],Momentum[p]]->p^2;

SEC2=SEC/.TBI[4,p^2,{{1,mz},{1,MChi}}]->b Bzc[s]\
/.TBI[4,p^2,{{1,MChi},{1,ma}}]->b Bac[s]\
/.TBI[4,p^2,{{1,mw},{1,MChi}}]->b Bwc[s]\
/.TAI[4,0,{{1,mw}}]->Aw\
/.TAI[4,0,{{1,MChi}}]->Ac\
/.TAI[4,0,{{1,ma}}]->Aa\
/.TAI[4,0,{{1,mz}}]->Az;
SEC2=Simplify[SEC2\
/. (Aw) -> 2 a p^2 B1cw[s] -(p^2 + MChi^2 - mw^2)*b Bwc[s] + Ac\
/. (Az) -> 2 a p^2 B1cz[s] -(p^2 + MChi^2 - mz^2)*b Bzc[s] + Ac\
/. (Aa) -> 2 a p^2 B1ca[s] -(p^2 + MChi^2 - ma^2)*b Bac[s] + Ac]




SigmaKN=-(e^2/(16*Pi^2*sw^2))*(4*c1*B1cw[s]+2);
SigmaMN=-(e^2*MChi/(16*Pi^2*sw^2))*(8*c0*Bwc[s]-4);
SENIbe = SigmaKN*p+SigmaMN
SigmaKC = -(e^2/(8*Pi^2 * sw^2)) * ( sw^2 *c1*B1ca[s] + cw^2 *c1*B1cz[s] + c1*B1cw[s] + 1 );
SigmaMC = -(e^2 MChi /(4*Pi^2 * sw^2)) * ( sw^2 *c0*Bac[s] + cw^2 *c0*Bzc[s] + c0*Bwc[s] - 1 );
SECIbe = (SigmaKC*p+SigmaMC)


(*a=-I
b=I*)
a=1;b=1;
Simplify[SENIbe-SEN2]


c1=I;
c0=-I;
Simplify[SECIbe-SEC2]

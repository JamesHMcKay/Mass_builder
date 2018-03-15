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


path = "/Users/jamesmckay/Documents/Programs/Mass_builder/";
kappa=1/(16\[Pi]^2);
(*STW=Sin[\[Theta]];
CTW=Cos[\[Theta]];*)

(* V5 *)

V5SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/VDM/output/math_data_V5_{1,2,3,4,5,6,7,8}_1.mx"}]];
V5SE = V5SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/VDM/output/math_data_V5_{9,10,11,12}_1.mx"}]];
V5SE = V5SE + SelfEnergyFinite*kappa;

SEc = makeFiniteAmplitude[V5SE,0,D]/.MVp->MChi/.MV0->MChi/.Pair[Momentum[p], Momentum[p]]->p^2;

SEc=SEc/.TBI[4,p^2,{{1,MZ},{1,MChi}}]-> I*Bzp\
/.TBI[4,p^2,{{1,MChi},{1,0}}]->I Bnp\
/.TBI[4,p^2,{{1,MChi},{1,MH}}]->I Bhp\
/.TBI[4,p^2,{{1,MW},{1,MChi}}]->I Bwp\
/.TAI[4,0,{{1,MW}}]-> I*Aw\
/.TAI[4,0,{{1,MChi}}]-> I*Ap\
/.TAI[4,0,{{1,MH}}]-> I*Ah\
/.TAI[4,0,{{1,MZ}}]-> I*Az;



(* V6 *)

V6SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/VDM/output/math_data_V6_{1,2,3,4,5}_1.mx"}]];
V6SE = V6SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/VDM/output/math_data_V6_{6,7,8}_1.mx"}]];
V6SE = V6SE + SelfEnergyFinite*kappa;

SEn = makeFiniteAmplitude[V6SE,0,D]/.MVp->MChi/.MV0->MChi/.Pair[Momentum[p], Momentum[p]]->p^2;
SEn=SEn/.TBI[4,p^2,{{1,MZ},{1,MChi}}]-> I*Bzp\
/.TBI[4,p^2,{{1,MChi},{1,0}}]->I Bnp\
/.TBI[4,p^2,{{1,MChi},{1,MH}}]->I Bhp\
/.TBI[4,p^2,{{1,MW},{1,MChi}}]->I Bwp\
/.TAI[4,0,{{1,MW}}]-> I*Aw\
/.TAI[4,0,{{1,MChi}}]-> I*Ap\
/.TAI[4,0,{{1,MH}}]-> I*Ah\
/.TAI[4,0,{{1,MZ}}]-> I*Az;


DeltaM = FullSimplify[-(SEc-SEn)/(2*MChi)/.p^2->MChi^2/.p^4->MChi^4,{CW^2+SW^2==1,CTW==MW/MZ}]


CBwp=FullSimplify[Coefficient[DeltaM,Bwp,1]]
CBzp=FullSimplify[Coefficient[DeltaM,Bzp,1]]
CBnp=FullSimplify[Coefficient[DeltaM,Bnp,1]]
CAz=FullSimplify[Coefficient[DeltaM,Az,1]]
CAw=FullSimplify[Coefficient[DeltaM,Aw,1]]
CAp=FullSimplify[Coefficient[DeltaM,Ap,1]]
FullSimplify[ DeltaM - (CBwp*Bwp+CBnp*Bnp+CBzp*Bzp + CAz*Az + CAw*Aw + CAp*Ap)]


F[m_] = -(30*MChi^4 + 26*MChi^2*m^2-5*m^4);
G[m_] = (12*MChi^2-5*m^2);
DeltaMTest = (EE^2 / (192*Pi^2*MChi^3*SW^2)) * ( Bwp*F[MW]-CW^2*Bzp*F[MZ] + Aw*G[MW]-CW^2*Az*G[MZ] 
+ 5*(MW^2-CW^2*MZ^2)*(Ap-2*MChi^2) + 30*SW^2*MChi^4*Bnp)


CW=Cos[\[Theta]];
SW=Sin[\[Theta]];
FullSimplify[DeltaM-DeltaMTest/.p->MChi,{CW^2+SW^2==1,CTW==MW/MZ}]


FullSimplify[SEn]


FullSimplify[Coefficient[SEn,p,2]]


FullSimplify[Coefficient[SEn,p,0]]


(* check counter-terms *)


V6SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/VDM/output/math_data_V6_{1,2,3,4,5}_1.mx"}]];
V6SE = V6SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/VDM/output/math_data_V6_{6,7,8}_1.mx"}]];
V6SE = V6SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/VDM/output/math_data_V6_1_1c.mx"}]];
V6SE = V6SE + SelfEnergyFinite*kappa;

SEn = makeFiniteAmplitude[V6SE,-1,D]/.MVp->MChi/.MV0->MChi/.MassBuilderCTM1->0/.MassBuilderCTZ1->0/.Pair[Momentum[p], Momentum[p]]->p^2;

eq1 = FullSimplify[Coefficient[SEn,p,2]];
eq2 = FullSimplify[Coefficient[SEn,p,0]];
solV6 = Solve[{eq1==0,eq2==0},{d6Z, d6M}];
Set @@@ solV6[[1]];


V5SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/VDM/output/math_data_V5_{1,2,3,4,5,6,7,8}_1.mx"}]];
V5SE = V5SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/VDM/output/math_data_V5_{9,10,11,12}_1.mx"}]];
V5SE = V5SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/VDM/output/math_data_V5_1_1c.mx"}]];
V5SE = V5SE + SelfEnergyFinite*kappa;

SEn = makeFiniteAmplitude[V5SE,-1,D]/.MVp->MChi/.MV0->MChi/.MassBuilderCTM1->0/.MassBuilderCTZ1->0/.Pair[Momentum[p], Momentum[p]]->p^2;

eq1 = FullSimplify[Coefficient[SEn,p,2]];
eq2 = FullSimplify[Coefficient[SEn,p,0]];
solV5 = Solve[{eq1==0,eq2==0},{d5Z, d5M}];
Set @@@ solV5[[1]];


FullSimplify[d6Z,{CW^2+SW^2==1}]
FullSimplify[d6M,{CW^2+SW^2==1}]
FullSimplify[d5Z,{CW^2+SW^2==1}]
FullSimplify[d5M,{CW^2+SW^2==1}]




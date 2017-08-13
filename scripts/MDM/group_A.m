(* ::Package:: *)

Quit[]


(* ::Section:: *)
(*Determine counter-term couplings for diagrams in group A*)


(* ::Subsection:: *)
(* Load required packages *)


$LoadTARCER = True;
$LoadFeynArts = True;
<< FeynCalc/FeynCalc.m;
AppendTo[$Path, "/Users/jamesmckay/Documents/Programs/Mass_builder/src/"];
<<MassBuilder.m;


SetOptions[DiracSlash,Dimension->D,FeynCalcInternal->True];SetOptions[DiracTrace,DiracTraceEvaluate->True];null=0;MassBuilderA[mass_,D_]:=TAI[D,0,{{1,mass}}];MassBuilderB[mass1_,mass2_,D_]:=TBI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2}}];MassBuilderJ[mass1_,mass2_,mass3_,D_]:=TJI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2},{1,mass3}}];MassBuilderT[mass1_,mass2_,mass3_,D_]:=TJI[D,Pair[Momentum[p,D],Momentum[p,D]],{{2,mass1},{1,mass2},{1,mass3}}];MassBuilderK[mass1_,mass2_,mass3_,D_]:=TJI[D,0,{{1,mass1},{1,mass2},{1,mass3}}];MassBuilderV[mass1_,mass2_,mass3_,mass4_,D_]:=TVI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2},{1,mass3},{1,mass4}}];MassBuilderF[mass1_,mass2_,mass3_,mass4_,mass5_,D_]:=TFI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2},{1,mass3},{1,mass4},{1,mass5}}];FCGV["EL"]=e;FCGV["SW"]=sw;FCGV["CW"]=cw;FCGV["MW"]=mw;FCGV["MZ"]=mz;FCGV["ME"]=me;FCGV["MM"]=mm;FCGV["ML"]=ml;FCGV["MU"]=mu;FCGV["MT"]=mt;FCGV["MD"]=md;FCGV["MS"]=ms;FCGV["MB"]=mb;ZNeu[1,1]=0;ZNeu[1,2]=1;ZNeu[1,3]=0;ZNeu[1,4]=0;ZNeu[2,1]=0;ZNeu[2,2]=0;ZNeu[2,3]=0;ZNeu[2,4]=0;ZNeu[3,1]=0;ZNeu[3,2]=0;ZNeu[3,3]=0;ZNeu[3,4]=0;ZNeu[4,1]=0;ZNeu[4,2]=0;ZNeu[4,3]=0;ZNeu[4,4]=0;UCha[1,1]=1;UCha[1,2]=0;UCha[2,1]=0;UCha[2,2]=0;VCha[1,1]=1;VCha[1,2]=0;VCha[2,1]=0;VCha[2,2]=0;CKM=IndexDelta;MCha[1]=MChi;MCha[2]=MChi;MNeu[1]=MChi;MNeu[2]=MChi;Mh0=mh;SBA=1;
SPD[p,p]=MChi^2;Pair[Momentum[p],Momentum[p]]=MChi^2;


(* ::Subsection:: *)
(* Get the amplitudes that have been computed previously*)
(* Change the path below to the location of your Mass_Builder/models directory *)


path = "/Users/jamesmckay/Documents/Programs/Mass_builder/";
kappa=1/(16\[Pi]^2);

(* Particle 1 *)
Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_23_2.mx"}]];
C23 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_24_2.mx"}]];
C24 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_25_2.mx"}]];
C25 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_26_2.mx"}]];
C26 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_27_2.mx"}]];
C27 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_28_2.mx"}]];
C28 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_29_2.mx"}]];
C29 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_30_2.mx"}]];
C30 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_31_2.mx"}]];
C31 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_32_2.mx"}]];
C32 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_33_2.mx"}]];
C33 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_34_2.mx"}]];
C34 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_35_2.mx"}]];
C35 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_36_2.mx"}]];
C36 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_37_2.mx"}]];
C37 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_38_2.mx"}]];
C38 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_39_2.mx"}]];
C39 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_40_2.mx"}]];
C40 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_41_2.mx"}]];
C41 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_42_2.mx"}]];
C42 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_43_2.mx"}]];
C43 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_11_2c.mx"}]];
C11c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_12_2c.mx"}]];
C12c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_13_2c.mx"}]];
C13c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_14_2c.mx"}]];
C14c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_15_2c.mx"}]];
C15c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_16_2c.mx"}]];
C16c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_17_2c.mx"}]];
C17c = SelfEnergyFinite*kappa^2;
Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_18_2c.mx"}]];
C18c = SelfEnergyFinite*kappa^2;

Csum1=makeFiniteAmplitude[C23+C24+C25+C26+C27+C28+C29+C30+C31+C32+C33+C34+C35+C36+C37+C38+C39+C40+C41+C42+C43,-1,D];
Csum2=makeFiniteAmplitude[C23+C24+C25+C26+C27+C28+C29+C30+C31+C32+C33+C34+C35+C36+C37+C38+C39+C40+C41+C42+C43,-2,D];

Cct1 = makeFiniteAmplitude[C11c+C12c+C13c+C14c+C15c+C16c+C17c+C18c,-1,D];
Cct2 = makeFiniteAmplitude[C11c+C12c+C13c+C14c+C15c+C16c+C17c+C18c,-2,D];

(* Doubly charged particle *)

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_17_2.mx"}]];
CC17 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_18_2.mx"}]];
CC18 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_19_2.mx"}]];
CC19 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_20_2.mx"}]];
CC20 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_21_2.mx"}]];
CC21 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_22_2.mx"}]];
CC22 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_23_2.mx"}]];
CC23 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_24_2.mx"}]];
CC24 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_25_2.mx"}]];
CC25 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_26_2.mx"}]];
CC26 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_27_2.mx"}]];
CC27 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_28_2.mx"}]];
CC28 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_29_2.mx"}]];
CC29 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_30_2.mx"}]];
CC30 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_31_2.mx"}]];
CC31 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_9_2c.mx"}]];
CC9c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_10_2c.mx"}]];
CC10c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_11_2c.mx"}]];
CC11c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_12_2c.mx"}]];
CC12c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_13_2c.mx"}]];
CC13c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_14_2c.mx"}]];
CC14c = SelfEnergyFinite*kappa^2;
CCsum1=makeFiniteAmplitude[CC17+CC18+CC19+CC20+CC21+CC22+CC23+CC24+CC25+CC26+CC27+CC28+CC29+CC30+CC31,-1,D];
CCsum2=makeFiniteAmplitude[CC17+CC18+CC19+CC20+CC21+CC22+CC23+CC24+CC25+CC26+CC27+CC28+CC29+CC30+CC31,-2,D];

CCct1 = makeFiniteAmplitude[CC9c+CC10c+CC11c+CC12c+CC13c+CC14c,-1,D];
CCct2 = makeFiniteAmplitude[CC9c+CC10c+CC11c+CC12c+CC13c+CC14c,-2,D];

(* Neutral particle *)

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_13_2.mx"}]];
N13 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_14_2.mx"}]];
N14 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_15_2.mx"}]];
N15 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_16_2.mx"}]];
N16 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_17_2.mx"}]];
N17 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_18_2.mx"}]];
N18 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_19_2.mx"}]];
N19 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_20_2.mx"}]];
N20 = SelfEnergyFinite*kappa^2;
Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_5_2c.mx"}]];
N5c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_6_2c.mx"}]];
N6c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_7_2c.mx"}]];
N7c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_8_2c.mx"}]];
N8c = SelfEnergyFinite*kappa^2;
Nsum1=makeFiniteAmplitude[N13+N14+N15+N16+N17+N18+N19+N20,-1,D];
Nsum2=makeFiniteAmplitude[N13+N14+N15+N16+N17+N18+N19+N20,-2,D];

Nct1 = makeFiniteAmplitude[N5c+N6c+N7c+N8c,-1,D];
Nct2 = makeFiniteAmplitude[N5c+N6c+N7c+N8c,-2,D];

{Nsum1,Nsum2,Csum1,Csum2,Nct1,Nct2,Cct1,Cct2,CCsum1,CCsum2,CCct1,CCct2}={Nsum1,Nsum2,Csum1,Csum2,Nct1,Nct2,Cct1,Cct2,CCsum1,CCsum2,CCct1,CCct2}/. (DiracGamma[Momentum[p, D], D].DiracGamma[6] + DiracGamma[Momentum[p, D], D].DiracGamma[7]) -> p\
/. MassBuilderP^2 -> Pair[Momentum[p],Momentum[p]]/. MassBuilderP -> Momentum[p]/. MassBuilderEpsilon -> \[Epsilon]/. MassBuilderZeta -> \[Zeta]/. MassBuilderAe[0] -> 0\
/. D -> 4 - 2 \[Epsilon]/. Pair[p, p] -> p^2/. D -> 4 - 2 \[Epsilon]/. Momentum[p] -> p\
/. MassBuilderB[MChi, mz] -> TBI[4, MChi^2, {{1, MChi}, {1, mz}}]\
/. MassBuilderB[mz, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, mz}}]\
/. MassBuilderB[MChi, mw] -> TBI[4, MChi^2, {{1, MChi}, {1, mw}}]\
/. MassBuilderB[mw, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, mw}}]\
/. MassBuilderB[ma, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, ma}}]\
/. MassBuilderB[MChi, ma] -> TBI[4, MChi^2, {{1, MChi}, {1, ma}}]\
/.p->MChi/. ma->0/.STW->sw/.CTW->cw;
(*sw=Sin[\[Theta]];
cw=Cos[\[Theta]];*)


dg2=16*g2^3*sw;
dg4=16*g2^3*cw;

eq1=Simplify[Coefficient[Simplify[Csum1+Cct1,{sw^2+cw^2==1}],TBI[4, MChi^2, {{1, MChi}, {1, 0}}]]]
eq2=Simplify[Coefficient[Simplify[Csum1+Cct1,{sw^2+cw^2==1}],TBI[4, MChi^2, {{1, MChi}, {1, mz}}]]]
eq3=Simplify[Coefficient[Simplify[Csum1+Cct1,{sw^2+cw^2==1}],TAI[4, 0, {1, mw}]]]
eq4=Simplify[Coefficient[Simplify[Csum1+Cct1,{sw^2+cw^2==1}],TAI[4, 0, {1, MChi}]]]
eq5=Simplify[Coefficient[Simplify[Csum1+Cct1,{sw^2+cw^2==1}],TAI[4, 0, {1, mz}]]]
eq6=Simplify[Coefficient[Simplify[Csum1+Cct1,{sw^2+cw^2==1}],TBI[4, MChi^2, {{1, MChi}, {1, mw}}]]]


dg5=32*g2^3*sw;
dg6=32*g2^3*cw;

dg3=dg7;
dg7=-32*g2^3/Sqrt[2];
eq12=Simplify[Coefficient[Simplify[CCsum1+CCct1,{sw^2+cw^2==1}],TAI[4, 0, {1, mw}]]]


dg1=dg8;
dg8=-48*g2^3/Sqrt[3];
eq14=Simplify[Coefficient[Simplify[Nsum1+Nct1,{sw^2+cw^2==1}],TAI[4, 0, {1, MChi}]]]








FullSimplify[Csum1+Cct1,{sw^2+cw^2==1}]


FullSimplify[Nsum1+Nct1,{sw^2+cw^2==1}]


FullSimplify[CCsum1+CCct1/.mz->mw/cw,{sw^2+cw^2==1}]


CForm[dg1]


CForm[dg2]


CForm[dg3]


CForm[dg4]


CForm[dg5]


CForm[dg6]


CForm[dg7]


CForm[dg8]

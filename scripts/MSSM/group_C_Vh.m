(* ::Package:: *)

Quit[]


(* ::Section:: *)
(*Determine counter-term couplings for diagrams in group C*)


(* ::Subsection:: *)
(* Load required packages *)


$LoadTARCER = True;
$LoadFeynArts = True;
<< FeynCalc/FeynCalc.m;
AppendTo[$Path, "/Users/jamesmckay/Documents/Programs/Mass_builder/src/"];
<<MassBuilder.m;


SetOptions[DiracSlash,Dimension->D,FeynCalcInternal->True];SetOptions[DiracTrace,DiracTraceEvaluate->True];null=0;MassBuilderA[mass_,D_]:=TAI[D,0,{{1,mass}}];MassBuilderB[mass1_,mass2_,D_]:=TBI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2}}];MassBuilderJ[mass1_,mass2_,mass3_,D_]:=TJI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2},{1,mass3}}];MassBuilderT[mass1_,mass2_,mass3_,D_]:=TJI[D,Pair[Momentum[p,D],Momentum[p,D]],{{2,mass1},{1,mass2},{1,mass3}}];MassBuilderK[mass1_,mass2_,mass3_,D_]:=TJI[D,0,{{1,mass1},{1,mass2},{1,mass3}}];MassBuilderV[mass1_,mass2_,mass3_,mass4_,D_]:=TVI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2},{1,mass3},{1,mass4}}];MassBuilderF[mass1_,mass2_,mass3_,mass4_,mass5_,D_]:=TFI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2},{1,mass3},{1,mass4},{1,mass5}}];FCGV["EL"]=e;FCGV["SW"]=sw;FCGV["CW"]=cw;FCGV["MW"]=mw;FCGV["MZ"]=mz;FCGV["ME"]=me;FCGV["MM"]=mm;FCGV["ML"]=ml;FCGV["MU"]=mu;FCGV["MT"]=mt;FCGV["MD"]=md;FCGV["MS"]=ms;FCGV["MB"]=mb;ZNeu[1,1]=0;ZNeu[1,2]=1;ZNeu[1,3]=0;ZNeu[1,4]=0;ZNeu[2,1]=0;ZNeu[2,2]=0;ZNeu[2,3]=0;ZNeu[2,4]=0;ZNeu[3,1]=0;ZNeu[3,2]=0;ZNeu[3,3]=0;ZNeu[3,4]=0;ZNeu[4,1]=0;ZNeu[4,2]=0;ZNeu[4,3]=0;ZNeu[4,4]=0;UCha[1,1]=1;UCha[1,2]=0;UCha[2,1]=0;UCha[2,2]=0;VCha[1,1]=1;VCha[1,2]=0;VCha[2,1]=0;VCha[2,2]=0;CKM=IndexDelta;MCha[1]=MChi;MCha[2]=MChi;MNeu[1]=MChi;MNeu[2]=MChi;Mh0=mh;SBA=1;


(* ::Subsection:: *)
(* Get the amplitudes that have been computed previously*)
(* Change the path below to the location of your Mass_Builder/models directory *)


path = "/Users/jamesmckay/Documents/Programs/Mass_builder/";
kappa=1/(16\[Pi]^2);
(*sw=Sin[\[Theta]];
cw=Cos[\[Theta]];*)
V3SE = 0;
ClearScalarProducts[];
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_1_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_2_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_3_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_4_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_5_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_6_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_13_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_14_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_15_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_16_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_17_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_18_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_19_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_20_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_21_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_22_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_23_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_24_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;
V3SEdiv = makeFiniteAmplitude[V3SE,-1,D];
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V3_1_1c.mx"}]];
V3ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V3ct+V3SEdiv,Pair[Momentum[p],Momentum[p]]]];
eq2 = FullSimplify[Coefficient[V3ct+V3SEdiv,Pair[Momentum[p],Momentum[p]],0]];
solV3 = Solve[{eq1==0,eq2==0},{dZW1, dMWsq1}]
(* Z boson *) 

V2SE = 0;
ClearScalarProducts[];
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_1_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_2_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_3_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_4_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_17_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_18_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_19_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_20_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_21_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_22_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_23_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_24_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_25_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;


V2SEdiv = makeFiniteAmplitude[V2SE,-1,D];
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V2_1_1c.mx"}]];
V2ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V2ct+V2SEdiv,Pair[Momentum[p],Momentum[p]]]];
eq2 = FullSimplify[Coefficient[V2ct+V2SEdiv,Pair[Momentum[p],Momentum[p]],0]];
solV2 = Solve[{eq1==0,eq2==0},{dMZsq1,dZZZ1}]

(* photon *) 
V1SE = 0;
ClearScalarProducts[];
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_1_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_2_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_12_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_13_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_14_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_15_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_16_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_17_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_18_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;
V1SEdiv = makeFiniteAmplitude[V1SE,-1,D];
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_1_1c.mx"}]];
V1ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V1ct+V1SEdiv,Pair[Momentum[p],Momentum[p]]]];
(*solV1 = Solve[{eq1\[Equal]0,eq2\[Equal]0},{dZAA1,dMAsq1}]*)
solV1 = Solve[{eq1==0},{dZAA1}]

(* photon Z mixing *) 

V12SE = 0;
ClearScalarProducts[];
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_1_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_2_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_12_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_13_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_14_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_15_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_16_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_17_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_18_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;
V12SEdiv = makeFiniteAmplitude[V12SE,-1,D];
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_V2_1_1c.mx"}]];
V12ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V12ct+V12SEdiv,Pair[Momentum[p],Momentum[p]]]];
eq2 = FullSimplify[Coefficient[V12ct+V12SEdiv,Pair[Momentum[p],Momentum[p]],0]];
solV12 = Solve[{eq1==0,eq2==0},{dZZA1,dZAZ1}]


Set @@@ solV1[[1]];
Set @@@ solV2[[1]];
Set @@@ solV3[[1]];
Set @@@ solV12[[1]];


SPD[p,p]=MChi^2;Pair[Momentum[p],Momentum[p]]=MChi^2;
(* Particle 1 *)

SelfEnergyTotal = 0;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_{1,2}_2.mx"}]];
N1 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_{3,4,9,10}_2.mx"}]];
N2 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_{5,6,11,12}_2.mx"}]];
N3 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_{7,8}_2.mx"}]];
N4 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_{31,32}_2.mx"}]];
N5 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_{39,40}_2.mx"}]];
N6 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_{41,42,45,46,49,50}_2.mx"}]];
N7 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_{43,44,47,48,51,52,57,58}_2.mx"}]];
N8 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_{53,54,59,60}_2.mx"}]];
N9 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_{55,56}_2.mx"}]];
N10 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_{3,4}_2c.mx"}]];
N3c = SelfEnergyFinite*kappa^2;

Nsum1a = makeFiniteAmplitude[N1+N2+N3+N4+N5+N6+N7+N8+N9+N10,-1,D];

Nsum1 = Nsum1a;

Nsum2a = makeFiniteAmplitude[N1+N2+N3+N4+N5+N6+N7+N8+N9+N10,-2,D];

Nsum2 = Nsum2a;

Nct1 = makeFiniteAmplitude[N3c,-1,D];
Nct2 = makeFiniteAmplitude[N3c,-2,D];

(* Particle 2 *)

SelfEnergyTotal = 0;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{1,10,84,91,92,104,105,115}_2.mx"}]];
C1 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{2,3}_2.mx"}]];
C2 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{4,5,6,13}_2.mx"}]];
C3 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{7,8,9,15,16}_2.mx"}]];
C4 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{11,12}_2.mx"}]];
C5 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_14_2.mx"}]];
C6 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_35_2.mx"}]];
C7 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{54,55}_2.mx"}]];
C8 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_68_2.mx"}]];
C9 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_75_2.mx"}]];
C10 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{85,86}_2.mx"}]];
C11 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_87_2.mx"}]];
C12 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{88,118}_2.mx"}]];
C13 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_89_2.mx"}]];
C14 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_90_2.mx"}]];
C15 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{93,94}_2.mx"}]];
C16 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{95,96}_2.mx"}]];
C35 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{97,98}_2.mx"}]];
C54 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{99,101,112,119}_2.mx"}]];
C55 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{100,113,120}_2.mx"}]];
C68 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "models/MSSM/output/math_data_F12_g1_102_2.mx"}]];
C75 = SelfEnergyFinite*kappa^2;
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_103_2.mx"}]];
C84 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "models/MSSM/output/math_data_F12_g1_{106,107,108,109}_2.mx"}]];
C85 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{110,111}_2.mx"}]];
C86 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_114_2.mx"}]];
C87 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{116,117}_2.mx"}]];
C88 = SelfEnergyFinite*kappa^2;


Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{4,5,6,7,8}_2c.mx"}]];
C8c = SelfEnergyFinite*kappa^2;
Csum1a = makeFiniteAmplitude[C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11+C12+C13+C14+C15+C16,-1,D];

Csum1b = makeFiniteAmplitude[C35+C54+C55+C68+C75+C84+C85+C86+C87+C88,-1,D];

Csum1 = Csum1a+Csum1b;
Csum2a = makeFiniteAmplitude[C1+C2+C3+C4+C5+C6+C7+C8+C9+C10+C11+C12+C13+C14+C15+C16,-2,D];
Csum2b = makeFiniteAmplitude[C35+C54+C55+C68+C75+C84+C85+C86+C87+C88,-2,D];

Csum2 = Csum2a+Csum2b;
Cct1 = makeFiniteAmplitude[C8c,-1,D];
Cct2 = makeFiniteAmplitude[C8c,-2,D];

{Csum1,Csum2,Nsum1,Nsum2,Cct1,Cct2,Nct1,Nct2}={Csum1,Csum2,Nsum1,Nsum2,Cct1,Cct2,Nct1,Nct2}/. (DiracGamma[Momentum[p, D], D].DiracGamma[6] + DiracGamma[Momentum[p, D], D].DiracGamma[7]) -> p\
/. MassBuilderP^2 -> Pair[Momentum[p],Momentum[p]]/. MassBuilderP -> Momentum[p]/. MassBuilderEpsilon -> \[Epsilon]/. MassBuilderZeta -> \[Zeta]/. MassBuilderAe[0] -> 0\
/. D -> 4 - 2 \[Epsilon]/. Pair[p, p] -> p^2/. D -> 4 - 2 \[Epsilon]/. Momentum[p] -> p\
/. MassBuilderB[MChi, mz] -> TBI[4, MChi^2, {{1, MChi}, {1, mz}}]\
/. MassBuilderB[mz, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, mz}}]\
/. MassBuilderB[MChi, mw] -> TBI[4, MChi^2, {{1, MChi}, {1, mw}}]\
/. MassBuilderB[mw, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, mw}}]\
/. MassBuilderB[ma, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, ma}}]\
/. MassBuilderB[MChi, ma] -> TBI[4, MChi^2, {{1, MChi}, {1, ma}}]\
/.p->MChi/.STW->sw/.CTW->cw;
sw=Sin[\[Theta]];
cw=Cos[\[Theta]];


Print["Difference between 1/epsilon^2 terms is: = "]
FullSimplify[Csum2+Cct2 -Nsum2-Nct2/.mz->mw/cw]


Print["Difference between 1/epsilon^2 terms is: = "]
Simplify[Csum1+Cct1 -Nsum1-Nct1 /.mz->mw/cw]


Simplify[Nsum1+Nct1]


CForm[FullSimplify[dZW1]]


CForm[FullSimplify[dMWsq1]]


CForm[FullSimplify[dMZsq1]]


CForm[FullSimplify[dZAA1]]


CForm[FullSimplify[dZZA1]]


(* ::InheritFromParent:: *)
(**)
(**)


CForm[FullSimplify[dZAZ1]]


CForm[FullSimplify[dZZZ1]]


(* Check that counter-term couplings agree with Ibe et al. (2013) *)
g=e/sw;
gammaZ = -e^2(-5/3);
ZgammaZ = (-e*g/cw)*(11/6-5*sw^2/3);
ZZ= (-g^2/cw^2)( -11/6+11*sw^2/3-5 sw^4/3);
WZ=-g^2 * (-11/6);

ZgammaM = (-e*g/cw)*(2-2 sw^2)*mz^2;
ZM = (-g^2/cw^2)*(-1+6 sw^2- 4 sw^4)*mz^2;
WM = (-g^2)*(-1+2 sw^2)*mz^2;


Print["Difference between counter-terms computed and those of Ibe et al."]
FullSimplify[gammaZ- dZAA1]

FullSimplify[ZZ-dZZZ1/.mz->mw/cw]
FullSimplify[WZ-dZW1/.mz->mw/cw]
FullSimplify[ZgammaZ+(1/2)(dZZA1+dZAZ1)/.mz->mw/cw]

FullSimplify[dMZsq1+dZZZ1*mz^2 + ZM /.mz->mw/cw]

FullSimplify[dMWsq1+dZW1*mw^2 + WM /.mz->mw/cw]

FullSimplify[ZgammaM-(1/2)dZZA1*mz^2/.mz->mw/cw]











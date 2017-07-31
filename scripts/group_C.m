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
Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_1_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_2_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_3_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_4_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_5_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_6_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_7_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_8_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_9_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_10_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_11_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_12_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_13_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_14_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_15_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_16_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_17_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_18_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

V3SEdiv = makeFiniteAmplitude[V3SE,-1,D];
Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V3_1_1c.mx"}]];
V3ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V3ct+V3SEdiv,Pair[Momentum[p],Momentum[p]]]];
eq2 = FullSimplify[Coefficient[V3ct+V3SEdiv,Pair[Momentum[p],Momentum[p]],0]];
solV3 = Solve[{eq1==0,eq2==0},{dZW1, dMWsq1}]
(* Z boson *) 

V2SE = 0;
ClearScalarProducts[];
Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_1_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_2_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_3_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_4_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_5_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_6_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_7_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_8_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_9_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_10_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_11_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_12_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_13_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_14_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

V2SEdiv = makeFiniteAmplitude[V2SE,-1,D];
Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V2_1_1c.mx"}]];
V2ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V2ct+V2SEdiv,Pair[Momentum[p],Momentum[p]]]];
eq2 = FullSimplify[Coefficient[V2ct+V2SEdiv,Pair[Momentum[p],Momentum[p]],0]];
solV2 = Solve[{eq1==0,eq2==0},{dMZsq1,dZZZ1}]

(* photon *) 
V1SE = 0;
ClearScalarProducts[];
Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_1_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_2_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_3_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_4_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_5_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_6_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_7_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_8_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_9_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

V1SEdiv = makeFiniteAmplitude[V1SE,-1,D];
Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_1_1c.mx"}]];
V1ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V1ct+V1SEdiv,Pair[Momentum[p],Momentum[p]]]];
(*solV1 = Solve[{eq1\[Equal]0,eq2\[Equal]0},{dZAA1,dMAsq1}]*)
solV1 = Solve[{eq1==0},{dZAA1}]

(* photon Z mixing *) 

V12SE = 0;
ClearScalarProducts[];
Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_V2_1_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_V2_2_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_V2_3_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_V2_4_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_V2_5_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_V2_6_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_V2_7_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_V2_8_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_V2_9_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

V12SEdiv = makeFiniteAmplitude[V12SE,-1,D];
Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_V1_V2_1_1c.mx"}]];
V12ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V12ct+V12SEdiv,Pair[Momentum[p],Momentum[p]]]];
eq2 = FullSimplify[Coefficient[V12ct+V12SEdiv,Pair[Momentum[p],Momentum[p]],0]];
solV12 = Solve[{eq1==0,eq2==0},{dZZA1,dZAZ1}]


Set @@@ solV1[[1]];
Set @@@ solV2[[1]];
Set @@@ solV3[[1]];
Set @@@ solV12[[1]];


Clear[dZZA1,dZAZ1,dZZZ1,dZAA1]


SPD[p,p]=MChi^2;Pair[Momentum[p],Momentum[p]]=MChi^2;
(* Particle 1 *)

SelfEnergyTotal = 0;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_1_2.mx"}]];
N1 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_2_2.mx"}]];
N2 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_3_2.mx"}]];
N3 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_4_2.mx"}]];
N4 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_5_2.mx"}]];
N5 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_6_2.mx"}]];
N6 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_7_2.mx"}]];
N7 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_8_2.mx"}]];
N8 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_9_2.mx"}]];
N9 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_10_2.mx"}]];
N10 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_11_2.mx"}]];
N11 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_12_2.mx"}]];
N12 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_13_2.mx"}]];
N13 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_14_2.mx"}]];
N14 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_15_2.mx"}]];
N15 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_16_2.mx"}]];
N16 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_32_2.mx"}]];
N32 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_33_2.mx"}]];
N33 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_34_2.mx"}]];
N34 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_35_2.mx"}]];
N35 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_45_2.mx"}]];
N45 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_46_2.mx"}]];
N46 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_47_2.mx"}]];
N47 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_48_2.mx"}]];
N48 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_49_2.mx"}]];
N49 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_50_2.mx"}]];
N50 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_51_2.mx"}]];
N51 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_52_2.mx"}]];
N52 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_53_2.mx"}]];
N53 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_54_2.mx"}]];
N54 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_55_2.mx"}]];
N55 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_56_2.mx"}]];
N56 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_57_2.mx"}]];
N57 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_58_2.mx"}]];
N58 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_59_2.mx"}]];
N59 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_60_2.mx"}]];
N60 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_61_2.mx"}]];
N61 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_62_2.mx"}]];
N62 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_63_2.mx"}]];
N63 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_64_2.mx"}]];
N64 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_65_2.mx"}]];
N65 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_66_2.mx"}]];
N66 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_67_2.mx"}]];
N67 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_68_2.mx"}]];
N68 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_69_2.mx"}]];
N69 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_70_2.mx"}]];
N70 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_71_2.mx"}]];
N71 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_72_2.mx"}]];
N72 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_73_2.mx"}]];
N73 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_74_2.mx"}]];
N74 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_4_2c.mx"}]];
N4c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_5_2c.mx"}]];
N5c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_6_2c.mx"}]];
N6c = SelfEnergyFinite*kappa^2;

NsumZZ1=makeFiniteAmplitude[N1+N2+N3+N10+N32+N33+N45+N46+N52+N51+N61+N62+N63+N70,-1,D];
NsumZZ2=makeFiniteAmplitude[N1+N2+N3+N10+N32+N33+N45+N46+N52+N51+N61+N62+N63+N70,-2,D];

NsumWW1=makeFiniteAmplitude[N4+N5+N6+N7+N8+N9+N11+N12+N13+N14+N15+N16+N34+N35+N47+N48+N49+N50+N53+N54+N55+N56+N57+N58+N59+N60+N64+N65+N66+N67+N68+N69+N71+N72+N73+N74,-1,D];
NsumWW2=makeFiniteAmplitude[N4+N5+N6+N7+N8+N9+N11+N12+N13+N14+N15+N16+N34+N35+N47+N48+N49+N50+N53+N54+N55+N56+N57+N58+N59+N60+N64+N65+N66+N67+N68+N69+N71+N72+N73+N74,-2,D];

Nct1 = makeFiniteAmplitude[N6c+N5c+N4c,-1,D];
Nct2 = makeFiniteAmplitude[N6c+N5c+N4c,-2,D];

(* Particle 2 *)

SelfEnergyTotal = 0;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_1_2.mx"}]];
C1 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_2_2.mx"}]];
C2 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_3_2.mx"}]];
C3 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_4_2.mx"}]];
C4 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_5_2.mx"}]];
C5 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_6_2.mx"}]];
C6 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_7_2.mx"}]];
C7 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_8_2.mx"}]];
C8 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_9_2.mx"}]];
C9 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_10_2.mx"}]];
C10 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_11_2.mx"}]];
C11 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_12_2.mx"}]];
C12 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_13_2.mx"}]];
C13 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_14_2.mx"}]];
C14 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_15_2.mx"}]];
C15 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_16_2.mx"}]];
C16 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_29_2.mx"}]];
C29 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_30_2.mx"}]];
C30 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_31_2.mx"}]];
C31 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_32_2.mx"}]];
C32 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_33_2.mx"}]];
C33 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_34_2.mx"}]];
C34 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "models/massiveMSSM/output/math_data_F12_g1_44_2.mx"}]];
C44 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_45_2.mx"}]];
C45 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_46_2.mx"}]];
C46 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_47_2.mx"}]];
C47 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_48_2.mx"}]];
C48 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_49_2.mx"}]];
C49 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_50_2.mx"}]];
C50 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "models/massiveMSSM/output/math_data_F12_g1_51_2.mx"}]];
C51 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "models/massiveMSSM/output/math_data_F12_g1_52_2.mx"}]];
C52 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_53_2.mx"}]];
C53 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_54_2.mx"}]];
C54 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_55_2.mx"}]];
C55 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_56_2.mx"}]];
C56 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_57_2.mx"}]];
C57 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_58_2.mx"}]];
C58 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_59_2.mx"}]];
C59 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_60_2.mx"}]];
C60 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_61_2.mx"}]];
C61 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_62_2.mx"}]];
C62 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_63_2.mx"}]];
C63 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "models/massiveMSSM/output/math_data_F12_g1_64_2.mx"}]];
C64 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "models/massiveMSSM/output/math_data_F12_g1_65_2.mx"}]];
C65 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_66_2.mx"}]];
C66 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_67_2.mx"}]];
C67 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_68_2.mx"}]];
C68 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_69_2.mx"}]];
C69 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_70_2.mx"}]];
C70 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_71_2.mx"}]];
C71 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_72_2.mx"}]];
C72 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_73_2.mx"}]];
C73 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_74_2.mx"}]];
C74 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "models/massiveMSSM/output/math_data_F12_g1_75_2.mx"}]];
C75 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_76_2.mx"}]];
C76 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_77_2.mx"}]];
C77 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_78_2.mx"}]];
C78 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_79_2.mx"}]];
C79 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_80_2.mx"}]];
C80 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_4_2c.mx"}]];
C4c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_5_2c.mx"}]];
C5c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_6_2c.mx"}]];
C6c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_7_2c.mx"}]];
C7c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_8_2c.mx"}]];
C8c = SelfEnergyFinite*kappa^2;

CsumZZ1=makeFiniteAmplitude[C4+C5+C6+C13+C32+C33+C47+C48+C57+C58+C63+C70+C71+C78,-1,D];
CsumZZ2=makeFiniteAmplitude[C4+C5+C6+C13+C32+C33+C47+C48+C57+C58+C63+C70+C71+C78,-2,D];

CsumWW1=makeFiniteAmplitude[C7+C8+C9+C14+C15+C16+C34+C49+C50+C59+C60+C61+C62+C72+C73+C74+C79+C80,-1,D];
CsumWW2=makeFiniteAmplitude[C7+C8+C9+C14+C15+C16+C34+C49+C50+C59+C60+C61+C62+C72+C73+C74+C79+C80,-2,D];

CsumZA1=makeFiniteAmplitude[C2+C3+C11+C12+C31+C30+C45+C46+C53+C54+C55+C56+C66+C67+C68+C69+C76+C77,-1,D];
CsumZA2=makeFiniteAmplitude[C2+C3+C11+C12+C31+C30+C45+C46+C53+C54+C55+C56+C66+C67+C68+C69+C76+C77,-2,D];

CsumAA1=makeFiniteAmplitude[C1+C10+C29+C44+C51+C52+C64+C65+C75,-1,D];
CsumAA2=makeFiniteAmplitude[C1+C10+C29+C44+C51+C52+C64+C65+C75,-2,D];

Cct1 = makeFiniteAmplitude[C4c+C5c+C6c+C7c+C8c,-1,D];
Cct2 = makeFiniteAmplitude[C4c+C5c+C6c+C7c+C8c,-2,D];


{CsumZZ1,CsumZZ2,CsumWW1,CsumWW2,NsumZZ1,NsumZZ2,NsumWW1,NsumWW2,CsumZA1,CsumZA2,CsumAA1,CsumAA2,Cct1,Cct2,Nct1,Nct2}={CsumZZ1,CsumZZ2,CsumWW1,CsumWW2,NsumZZ1,NsumZZ2,NsumWW1,NsumWW2,CsumZA1,CsumZA2,CsumAA1,CsumAA2,Cct1,Cct2,Nct1,Nct2}/. (DiracGamma[Momentum[p, D], D].DiracGamma[6] + DiracGamma[Momentum[p, D], D].DiracGamma[7]) -> p\
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


FullSimplify[CsumZZ2+CsumWW2+CsumZA2+CsumAA2+Cct2 -NsumZZ2-NsumWW2-Nct2/.mz->mw/cw]


Simplify[CsumZZ1+CsumWW1+CsumZA1+CsumAA1+Cct1-(NsumWW1+NsumZZ1+Nct1) /.mz->mw/cw]


CForm[FullSimplify[dZW1]]


CForm[FullSimplify[dMWsq1]]


CForm[FullSimplify[dMZsq1]]


CForm[FullSimplify[dZAA1]]


CForm[FullSimplify[dZZA1]]


(* ::InheritFromParent:: *)
(**)


CForm[FullSimplify[dZAZ1]]


CForm[FullSimplify[dZZZ1]]


FullSimplify[dZW1]


(* Check that counter-term couplings agree with Ibe et al. (2013) *)
g=e/sw;
gammaZ = -e^2(-5/3);
ZgammaZ = (-e*g/cw)*(11/6-5*sw^2/3);
ZZ= (-g^2/cw^2)( -11/6+11*sw^2/3-5 sw^4/3);
WZ=-g^2 * (-11/6);

ZgammaM = (-e*g/cw)*(2-2 sw^2)*mz^2;
ZM = (-g^2/cw^2)*(-1+6 sw^2- 4 sw^4)*mz^2;
WM = (-g^2)*(-1+2 sw^2)*mz^2;


FullSimplify[gammaZ- dZAA1]

FullSimplify[ZZ-dZZZ1/.mz->mw/cw]
FullSimplify[WZ-dZW1/.mz->mw/cw]
FullSimplify[ZgammaZ+(1/2)(dZZA1+dZAZ1)/.mz->mw/cw]

FullSimplify[dMZsq1+dZZZ1*mz^2 + ZM /.mz->mw/cw]

FullSimplify[dMWsq1+dZW1*mw^2 + WM /.mz->mw/cw]

FullSimplify[ZgammaM-(1/2)dZZA1*mz^2/.mz->mw/cw]





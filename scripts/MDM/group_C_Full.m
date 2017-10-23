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
STW=Sin[\[Theta]];
S2TW=Sin[2*\[Theta]];
CTW=Cos[\[Theta]];
C2TW=Cos[2*\[Theta]];
V3SE = 0;
ClearScalarProducts[];
Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_1_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_2_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_3_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_4_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_5_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_6_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_7_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_8_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_9_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_10_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_11_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_12_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_13_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_14_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_15_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_16_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_17_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_18_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_19_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_20_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_21_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_22_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_23_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_24_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;


Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_25_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_26_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_27_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_28_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_29_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_30_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_31_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_32_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_33_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_34_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_35_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_36_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_37_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

V3SEdiv = makeFiniteAmplitude[V3SE,-1,D];
Get[FileNameJoin[{path, "/models/MDM/output/math_data_V3_1_1c.mx"}]];
V3ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V3ct+V3SEdiv,Pair[Momentum[p],Momentum[p]]]];
eq2 = FullSimplify[Coefficient[V3ct+V3SEdiv,Pair[Momentum[p],Momentum[p]],0]];
solV3 = Solve[{eq1==0,eq2==0},{dZW1, dMWsq1}]
(* Z boson *) 

V2SE = 0;
ClearScalarProducts[];
Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_1_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_2_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_3_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_4_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_5_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_6_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_7_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_8_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_9_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_10_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_11_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_12_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_13_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_14_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_15_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_16_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_17_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_18_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_19_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_20_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_21_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_22_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_23_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_24_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_25_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_26_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

V2SEdiv = makeFiniteAmplitude[V2SE,-1,D];
Get[FileNameJoin[{path, "/models/MDM/output/math_data_V2_1_1c.mx"}]];
V2ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V2ct+V2SEdiv,Pair[Momentum[p],Momentum[p]]]/.g1->g2*STW/CTW/.v->2*mw/g2/.mz->mw/CTW];
eq2 = FullSimplify[Coefficient[V2ct+V2SEdiv,Pair[Momentum[p],Momentum[p]],0]/.g1->g2*STW/CTW/.v->2*mw/g2/.mz->mw/CTW];
solV2 = Solve[{eq1==0,eq2==0},{dMZsq1,dZZZ1}]


(* photon *) 
V1SE = 0;
ClearScalarProducts[];
Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_1_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_2_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_3_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_4_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_5_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_6_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_7_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_8_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_9_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_10_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_11_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_12_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_13_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_14_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_15_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_16_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_17_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_18_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_19_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;


V1SEdiv = makeFiniteAmplitude[V1SE,-1,D];
Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_1_1c.mx"}]];
V1ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V1ct+V1SEdiv,Pair[Momentum[p],Momentum[p]]]/.g1->g2*STW/CTW];
(*solV1 = Solve[{eq1\[Equal]0,eq2\[Equal]0},{dZAA1,dMAsq1}]*)
solV1 = Solve[{eq1==0},{dZAA1}]



(* photon Z mixing *) 

V12SE = 0;
ClearScalarProducts[];
Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_1_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_2_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_3_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_4_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_5_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_6_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_7_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_8_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_9_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_10_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_11_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_12_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_13_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_14_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_15_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_16_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_17_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_18_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;


Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_19_1.mx"}]];
V12SE = V12SE + SelfEnergyFinite*kappa;


V12SEdiv = makeFiniteAmplitude[V12SE,-1,D];
Get[FileNameJoin[{path, "/models/MDM/output/math_data_V1_V2_1_1c.mx"}]];
V12ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[V12ct+V12SEdiv,Pair[Momentum[p],Momentum[p]]]/.g1->g2*STW/CTW/.v->2*mw/g2];
eq2 = FullSimplify[Coefficient[V12ct+V12SEdiv,Pair[Momentum[p],Momentum[p]],0]/.g1->g2*STW/CTW/.v->2*mw/g2/.mz->mw/CTW];
solV12 = Solve[{eq1==0,eq2==0},{dZZA1,dZAZ1}]



Set @@@ solV1[[1]];
Set @@@ solV2[[1]];
Set @@@ solV3[[1]];
Set @@@ solV12[[1]];


SPD[p,p]=MChi^2;Pair[Momentum[p],Momentum[p]]=MChi^2;
(* Particle 1 *)


SEtotal=0;
SEtotalc=0;
 
Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{1,2,3,14,15}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{4,5,6,16}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{7,8,9,10,11,12}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_13_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{17,18}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{19,20}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{21,22}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 

 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{44,45,46}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_47_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_48_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_49_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_50_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_51_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_52_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{53,54}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{55,56,57,58,59,60}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{61,64}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{62,65}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{63,66}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{67,70}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{68,71}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{69,72}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{73,74,75,76}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{77,78,79,80,81,82}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_83_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_84_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_85_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_86_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_87_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_88_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{89,90}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{91,92,93,94,95,96,97,98,99}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_100_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_101_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_102_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_103_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_104_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_105_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_106_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_107_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_108_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_109_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{110,111,112,113,114,115,116,117,118}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_119_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_120_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_121_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_122_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_123_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_124_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_125_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_126_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_127_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{128,129,130}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_144_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{145,146}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{147,148}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{149,150,151,152}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{153,154}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{155,156,157,158}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{159,160}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{161,162}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{163,164}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{165,166}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{167,168}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_169_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{170,171}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{172,173,174,175}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{176,177}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{178,179}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{180,181}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{182,183}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_184_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{185,186}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_187_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{188,189}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{190,191}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 
 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_5_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{6,7}_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_8_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{9,10}_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa;


Csum1=makeFiniteAmplitude[SEtotal,-1,D];
Csum2=makeFiniteAmplitude[SEtotal,-2,D];

Cct1 = makeFiniteAmplitude[SEtotalc,-1,D];
Cct2 = makeFiniteAmplitude[SEtotalc,-2,D];

SEtotal=0;
SEtotalc=0;

(* Doubly charged particle *)

Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{1,2,3}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{4,5,6}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{7,8,9}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{10,11,12}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_13_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{14,15,16}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 
 
 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{32,33,34}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_35_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_36_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_37_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_38_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_39_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_40_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{41,42}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{43,44,45,46,47,48}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{49,52}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{50,53}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{51,54}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{55,58}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{56,59}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{57,60}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{61,62,63,64}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{65,66,67,68,69,70}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_71_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_72_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_73_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_74_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_75_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_76_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{77,78}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{79,80,81,82,83,84,85,86,87}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_88_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_89_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_90_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_91_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_92_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_93_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 
  Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_94_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_95_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_96_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{97,98}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 
 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_109_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{110,111}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_112_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_113_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_114_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_115_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{116,117}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{118,119,120,121}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{122,123}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{124,126}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{125,127}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_128_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{129,130}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{131,132,133,134}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{135,136}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_137_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_138_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_139_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_140_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{141,142}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_143_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_144_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_145_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 


 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_4_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{5,6}_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_7_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_8_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 

CCsum1=makeFiniteAmplitude[SEtotal,-1,D];
CCsum2=makeFiniteAmplitude[SEtotal,-2,D];

CCct1 = makeFiniteAmplitude[SEtotalc,-1,D];
CCct2 = makeFiniteAmplitude[SEtotalc,-2,D];
(* Neutral particle *)

SEtotal=0;
SEtotalc=0;

Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{1,2}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{3,4}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{5,6,11,12}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{7,8}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{9,10}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 
 
 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{39,48}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{40,49}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{41,50}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{42,51}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{43,52}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{44,53}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{45,54}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{46,55}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{47,56}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{57,58,59,60}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{69,70,71,72}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{73,74}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{75,76}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{77,78}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{79,80}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{81,82}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{83,84}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{85,86}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{87,88}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{89,90}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 


 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{3,4}_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 

Nsum1=makeFiniteAmplitude[SEtotal,-1,D];
Nsum2=makeFiniteAmplitude[SEtotal,-2,D];

Nct1 = makeFiniteAmplitude[SEtotalc,-1,D];
Nct2 = makeFiniteAmplitude[SEtotalc,-2,D];




{Csum1,Csum2,Nsum1,Nsum2,Cct1,Cct2,Nct1,Nct2}={Csum1,Csum2,Nsum1,Nsum2,Cct1,Cct2,Nct1,Nct2}/. (DiracGamma[Momentum[p, D], D].DiracGamma[6] + DiracGamma[Momentum[p, D], D].DiracGamma[7]) -> p\
/. MassBuilderP^2 -> Pair[Momentum[p],Momentum[p]]/. MassBuilderP -> Momentum[p]/. MassBuilderEpsilon -> \[Epsilon]/. MassBuilderZeta -> \[Zeta]/. MassBuilderAe[0] -> 0\
/. D -> 4 - 2 \[Epsilon]/. Pair[p, p] -> p^2/. D -> 4 - 2 \[Epsilon]/. Momentum[p] -> p\
/. MassBuilderB[a_, b_] -> TBI[4, MChi^2, {{1, a}, {1, b}}]\
/.p->MChi/.STW->sw/.CTW->cw;
sw=Sin[\[Theta]];
cw=Cos[\[Theta]];


Simplify[Csum2+Cct2 -Nsum2-Nct2/.mz->mw/cw/.g1->g2*STW/CTW/.v->2*mw/g2/.mz->mw/CTW]


Simplify[ Nsum2+Nct2/.mz->mw/cw/.g1->g2*STW/CTW/.v->2*mw/g2/.mz->mw/CTW]


Simplify[ Coefficient[Nsum2+Nct2,mf,2]/.mz->mw/cw/.g1->g2*STW/CTW/.v->2*mw/g2/.mz->mw/CTW,TimeConstraint->100000]


Simplify[Csum1+Cct1/.mz->mw/cw,TimeConstraint->100000]


 Simplify[Coefficient[Csum1+Cct1/.mz->mw/cw,TBI[4, MChi^2, {{1, MChi}, {1, mz}}]]/.g1->g2*STW/CTW/.v->2*mw/g2/.mz->mw/CTW]


Simplify[Coefficient[Nsum1+Nct1/.mz->mw/cw,TBI[4, MChi^2, {{1, MChi}, {1, ma}}]]]


Csum2+Cct2


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





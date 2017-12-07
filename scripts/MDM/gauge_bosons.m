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


(* ::Subsection:: *)
(* Get the amplitudes that have been computed previously*)
(* Change the path below to the location of your Mass_Builder/models directory *)


path = "/Users/jamesmckay/Documents/Programs/Mass_builder/";
kappa=1/(16\[Pi]^2);
STW=Sin[\[Theta]];
sw=Sin[\[Theta]];
CTW=Cos[\[Theta]];
cw=Cos[\[Theta]];

SEtotal=0;
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_V1_12_1.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa;
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_V1_13_1.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa;
SEMDM=makeFiniteAmplitude[SEtotal,0,D];
SEtotal=0;
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MSSM/output/math_data_V1_12_1.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa;
SEMSSM=makeFiniteAmplitude[SEtotal,0,D];


FullSimplify[SEMSSM*I/.g2->e/sw]


FullSimplify[SEMDM*I/.g2->e/sw]


SEtotal=0;
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_V2_17_1.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa;
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_V2_18_1.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa;
SEMDM=makeFiniteAmplitude[SEtotal,0,D];
SEtotal=0;
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MSSM/output/math_data_V2_17_1.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa;
SEMSSM=makeFiniteAmplitude[SEtotal,0,D];


FullSimplify[SEMSSM*I/.g2->e/sw]


FullSimplify[SEMDM*I/.g2->e/sw]


SEtotal=0;
 (*Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_V3_25_1.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa;*)
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_V3_26_1.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa;
SEMDM=makeFiniteAmplitude[SEtotal,0,D];
SEtotal=0;
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MSSM/output/math_data_V3_13_1.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa;
SEMSSM=makeFiniteAmplitude[SEtotal,0,D];


FullSimplify[SEMSSM*I/.g2->e/sw]


FullSimplify[SEMDM*I/.g2->e/sw]


SEtotal=0;
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_V1_V2_12_1.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa;
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_V1_V2_13_1.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa;
SEMDM=makeFiniteAmplitude[SEtotal,0,D];
SEtotal=0;
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MSSM/output/math_data_V1_V2_12_1.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa;
SEMSSM=makeFiniteAmplitude[SEtotal,0,D];


FullSimplify[SEMSSM*I/.g2->e/sw]


FullSimplify[SEMDM*I/.g2->e/sw]


SEtotal=0;
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_V3_25_1.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa;

SEMDM=makeFiniteAmplitude[SEtotal,0,D];
SEtotal=0;
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_V3_26_1.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa;
SEMSSM=makeFiniteAmplitude[SEtotal,0,D];


FullSimplify[SEMSSM*I/.g2->e/sw/.mb->mt]





FullSimplify[SEMDM*I/.g2->e/sw]

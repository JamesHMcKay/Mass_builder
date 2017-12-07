(* ::Package:: *)

Quit[]


(* ::Section:: *)
(*Determine counter-term couplings for diagrams in group B*)


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
CTW=Cos[\[Theta]];

(* F5 *)

F5SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_{1,2,3,4}_1.mx"}]];
F5SE = F5SE + SelfEnergyFinite*kappa;

F5SE = makeFiniteAmplitude[F5SE,0,D];
F5SEdiv = makeFiniteAmplitude[F5SE,-1,D];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F5_1_1c.mx"}]];

F5ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;

eq1 = FullSimplify[Coefficient[F5ct+F5SEdiv,p]];
eq2 = FullSimplify[Coefficient[F5ct+F5SEdiv,p,0]];
solF5 = Solve[{eq1==0,eq2==0},{d5Z, d5m}]

(* F6 *)

F6SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_{1,2,3}_1.mx"}]];
F6SE = F6SE + SelfEnergyFinite*kappa;
F6SE = makeFiniteAmplitude[F6SE,0,D];

F6SEdiv = makeFiniteAmplitude[F6SE,-1,D];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F6_1_1c.mx"}]];

F6ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;

eq1 = FullSimplify[Coefficient[F6ct+F6SEdiv,p]];
eq2 = FullSimplify[Coefficient[F6ct+F6SEdiv,p,0]];
solF6 = Solve[{eq1==0,eq2==0},{d6Z, d6m}]
(* F7 *)

F7SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_{1,2}_1.mx"}]];
F7SE = F7SE + SelfEnergyFinite*kappa;

F7SEdiv = makeFiniteAmplitude[F7SE,-1,D];

F7SE = makeFiniteAmplitude[F7SE,0,D];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_F7_1_1c.mx"}]];

F7ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;

eq1 = FullSimplify[Coefficient[F7ct+F7SEdiv,p]];
eq2 = FullSimplify[Coefficient[F7ct+F7SEdiv,p,0]];
solF7 = Solve[{eq1==0,eq2==0},{d7Z, d7m}]
Set @@@ solF5[[1]];
Set @@@ solF6[[1]];
Set @@@ solF7[[1]];


SPD[p,p]=MChi^2;Pair[Momentum[p],Momentum[p]]=MChi^2;
(* Particle 1 *)
SelfEnergyTotal = 0;



F5SE-F6SE
F5SE-F7SE




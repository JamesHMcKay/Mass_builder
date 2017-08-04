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


Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_13_2.mx"}]];
N13 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_14_2.mx"}]];
N14 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_15_2.mx"}]];
N15 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_16_2.mx"}]];
N16 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_17_2.mx"}]];
N17 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_18_2.mx"}]];
N18 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_5_2c.mx"}]];
N5c =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_6_2c.mx"}]];
N6c =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_7_2c.mx"}]];
N7c =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_8_2c.mx"}]];
N8c =SelfEnergyFinite*kappa^2;

sumN1=makeFiniteAmplitude[N13+N14+N15+N16+N17+N18+N5c+N6c+N7c+N8c,-1,D];

sumN2=makeFiniteAmplitude[N13+N14+N15+N16+N17+N18+N5c+N6c+N7c+N8c,-2,D];

(* Particle 2 *)
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_17_2.mx"}]];
C17 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_18_2.mx"}]];
C18 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_19_2.mx"}]];
C19 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_20_2.mx"}]];
C20 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_21_2.mx"}]];
C21 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_22_2.mx"}]];
C22 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_23_2.mx"}]];
C23 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_24_2.mx"}]];
C24 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_25_2.mx"}]];
C25 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_9_2c.mx"}]];
C9c =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_10_2c.mx"}]];
C10c =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_11_2c.mx"}]];
C11c =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_12_2c.mx"}]];
C12c =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_13_2c.mx"}]];
C13c =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_14_2c.mx"}]];
C14c =SelfEnergyFinite*kappa^2;


sumC1=makeFiniteAmplitude[C17+C18+C19+C20+C21+C22+C23+C24+C25+C9c+C10c+C11c+C12c+C13c+C14c,-1,D];


sumC2=makeFiniteAmplitude[C17+C18+C19+C20+C21+C22+C23+C24+C25+C9c+C10c+C11c+C12c+C13c+C14c,-2,D];



{sumC1,sumN1,sumC2,sumN2}={sumC1,sumN1,sumC2,sumN2}/. (DiracGamma[Momentum[p, D], D].DiracGamma[6] + DiracGamma[Momentum[p, D], D].DiracGamma[7]) -> p\
/. MassBuilderP^2 -> Pair[Momentum[p],Momentum[p]]/. MassBuilderP -> Momentum[p]/. MassBuilderEpsilon -> \[Epsilon]/. MassBuilderZeta -> \[Zeta]/. MassBuilderAe[0] -> 0\
/. D -> 4 - 2 \[Epsilon]/. Pair[p, p] -> p^2/. D -> 4 - 2 \[Epsilon]/. Momentum[p] -> p\
/. MassBuilderB[MChi, mz] -> TBI[4, MChi^2, {{1, MChi}, {1, mz}}]\
/. MassBuilderB[mz, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, mz}}]\
/. MassBuilderB[MChi, mw] -> TBI[4, MChi^2, {{1, MChi}, {1, mw}}]\
/. MassBuilderB[mw, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, mw}}]\
/. MassBuilderB[ma, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, ma}}]\
/. MassBuilderB[MChi, ma] -> TBI[4, MChi^2, {{1, MChi}, {1, ma}}]\
/.p->MChi/. ma->0/.STW->sw/.CTW->cw;


eq6=Coefficient[FullSimplify[sumN1,{sw^2+cw^2==1,Cw1==-2,Cw2==-2}],TAI[4, 0, {1, mw}]]
sol1=Solve[eq6==0,{dg2}]


Set@@@sol1[[1]];


eq2=Coefficient[FullSimplify[sumC1,{sw^2+cw^2==1,Cw1+Cw2==-4}],TBI[4, MChi^2, {{1, MChi}, {1, 0}}]]
eq4=Coefficient[FullSimplify[sumC1,{sw^2+cw^2==1,Cw1+Cw2==-4}],TBI[4, MChi^2, {{1, MChi}, {1, mz}}]]
eq5=Coefficient[FullSimplify[sumC1,{sw^2+cw^2==1,Cw1+Cw2==-4}],TAI[4, 0, {1, MChi}]]
eq7=Coefficient[FullSimplify[sumC1,{sw^2+cw^2==1,Cw1+Cw2==-4}],TAI[4, 0, {1, mz}]]


sol2=Solve[{eq2==0,eq4==0,eq5==0,eq7==0},{Ca,Cz2}]


Set@@@sol2[[1]];

(*Cw1+Cw2==-4;*)

sw=Sin[Theta];
cw=Cos[Theta];


eq1=FullSimplify[sumN1-sumC1/.dg2->4*e^3/sw^3,{sw^2+cw^2==1,Cw1+Cw2==-4}]
eq1=FullSimplify[sumN2-sumC2/.dg2->4*e^3/sw^3,{sw^2+cw^2==1,Cw1+Cw2==-4}]


eq1=FullSimplify[sumN1/.dg2->4*e^3/sw^3,{sw^2+cw^2==1,Cw1+Cw2==-4}]
eq1=FullSimplify[sumN2/.dg2->4*e^3/sw^3,{sw^2+cw^2==1,Cw1+Cw2==-4}]
eq1=FullSimplify[sumC1/.dg2->4*e^3/sw^3,{sw^2+cw^2==1,Cw1+Cw2==-4}]
eq1=FullSimplify[sumC2/.dg2->4*e^3/sw^3,{sw^2+cw^2==1,Cw1+Cw2==-4}]


dg2*kappa

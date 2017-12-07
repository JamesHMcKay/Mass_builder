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
sw=Sin[\[Theta]];
cw=Cos[\[Theta]];

(* F11 *)

F11SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_{1,2}_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

F11SEdiv = makeFiniteAmplitude[F11SE,-1,D];

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_1_1c.mx"}]];

F11ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;

eq1 = FullSimplify[Coefficient[F11ct+F11SEdiv,p]];
eq2 = FullSimplify[Coefficient[F11ct+F11SEdiv,p,0]];
solF11 = Solve[{eq1==0,eq2==0},{d2Z, d2m}]



(* F12 *)

F12SE = 0;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_{1,2,3}_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

F12SEdiv = makeFiniteAmplitude[F12SE,-1,D];
Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_1_1c.mx"}]];
F12ct= makeFiniteAmplitude[SelfEnergyFinite*kappa,-1,D]/.MassBuilderCTM1->0/.MassBuilderCTZ1->0;
eq1 = FullSimplify[Coefficient[F12ct+F12SEdiv,p]];
eq2 = FullSimplify[Coefficient[F12ct+F12SEdiv,p,0]];
solF12 = Solve[{eq1==0,eq2==0},{d1Z, d1m}]




Set @@@ solF11[[1]];
Set @@@ solF12[[1]];


SPD[p,p]=MChi^2;Pair[Momentum[p],Momentum[p]]=MChi^2;
(* Particle 1 *)

SelfEnergyTotal = 0;



Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_33_2.mx"}]];
N33 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_34_2.mx"}]];
N34 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_35_2.mx"}]];
N35 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_36_2.mx"}]];
N36 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_37_2.mx"}]];
N37 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_38_2.mx"}]];
N38 = SelfEnergyFinite*kappa^2;





Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_1_2c.mx"}]];
N1c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F11_g1_2_2c.mx"}]];
N2c = SelfEnergyFinite*kappa^2;

Nsum1=makeFiniteAmplitude[N33+N34+N35+N36+N37+N38,-1,D];
Nsum2=makeFiniteAmplitude[N33+N34+N35+N36+N37+N38,-2,D];

Nct1 = makeFiniteAmplitude[N1c+N2c,-1,D];
Nct2 = makeFiniteAmplitude[N1c+N2c,-2,D];

(* Particle 2 *)

SelfEnergyTotal = 0;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_76_2.mx"}]];
C76 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_77_2.mx"}]];
C77 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_78_2.mx"}]];
C78 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_79_2.mx"}]];
C79 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_80_2.mx"}]];
C80 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_81_2.mx"}]];
C81 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_82_2.mx"}]];
C82 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_83_2.mx"}]];
C83 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_1_2c.mx"}]];
C1c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_2_2c.mx"}]];
C2c = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_F12_g1_3_2c.mx"}]];
C3c = SelfEnergyFinite*kappa^2;

Csum1=makeFiniteAmplitude[C76+C77+C78+C79+C80+C81+C82+C83,-1,D];
Csum2=makeFiniteAmplitude[C76+C77+C78+C79+C80+C81+C82+C83,-2,D];

Cct1 = makeFiniteAmplitude[C1c+C2c+C3c,-1,D];
Cct2 = makeFiniteAmplitude[C1c+C2c+C3c,-2,D];

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


FullSimplify[Nsum2+Nct2-(Csum2+Cct2)/.mz->mw/cw]


FullSimplify[Nsum1+Nct1-(Csum1+Cct1)/.mz->mw/cw]


CForm[FullSimplify[d1Z]]


CForm[FullSimplify[d1m]]


CForm[FullSimplify[d2Z]]


CForm[FullSimplify[d2m]]


FullSimplify[d1m*kappa]
FullSimplify[d2m*kappa]
FullSimplify[d1Z*kappa]
FullSimplify[d2Z*kappa]




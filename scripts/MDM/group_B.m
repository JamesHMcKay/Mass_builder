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



SEtotal=0;
SEtotalc=0;
 

 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_131_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{132,133}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_134_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{135,136,137}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{138,139,140}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{141,142,143}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 

 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{1,2,3,4}_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 

Csum1=makeFiniteAmplitude[SEtotal,-1,D];
Csum2=makeFiniteAmplitude[SEtotal,-2,D];

Cct1 = makeFiniteAmplitude[SEtotalc,-1,D];
Cct2 = makeFiniteAmplitude[SEtotalc,-2,D];

SEtotal=0;
SEtotalc=0;

(* Doubly charged particle *)

 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_99_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{100,101}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_102_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{103,104}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{105,106}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{107,108}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 

 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{1,2,3}_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 

CCsum1=makeFiniteAmplitude[SEtotal,-1,D];
CCsum2=makeFiniteAmplitude[SEtotal,-2,D];

CCct1 = makeFiniteAmplitude[SEtotalc,-1,D];
CCct2 = makeFiniteAmplitude[SEtotalc,-2,D];
(* Neutral particle *)

SEtotal=0;
SEtotalc=0;

 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{61,62}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{63,64}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{65,66,67,68}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 

 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{1,2}_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 


Nsum1=makeFiniteAmplitude[SEtotal,-1,D];
Nsum2=makeFiniteAmplitude[SEtotal,-2,D];

Nct1 = makeFiniteAmplitude[SEtotalc,-1,D];
Nct2 = makeFiniteAmplitude[SEtotalc,-2,D];


{Csum1,Csum2,Nsum1,Nsum2,Cct1,Cct2,Nct1,Nct2,CCsum1,CCsum2,CCct1,CCct2}={Csum1,Csum2,Nsum1,Nsum2,Cct1,Cct2,Nct1,Nct2,CCsum1,CCsum2,CCct1,CCct2}/. (DiracGamma[Momentum[p, D], D].DiracGamma[6] + DiracGamma[Momentum[p, D], D].DiracGamma[7]) -> p\
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
FullSimplify[Nsum2+Nct2-(CCsum2+CCct2)/.mz->mw/cw]


FullSimplify[Nsum2+Nct2/.mz->mw/cw]
FullSimplify[Csum2+Cct2/.mz->mw/cw]
FullSimplify[CCsum2+CCct2/.mz->mw/cw]


FullSimplify[Nsum1+Nct1-(Csum1+Cct1)/.mz->mw/cw]
FullSimplify[Nsum1+Nct1-(CCsum1+CCct1)/.mz->mw/cw]


CForm[FullSimplify[d5Z]]


CForm[FullSimplify[d5m]]


CForm[FullSimplify[d6Z]]


CForm[FullSimplify[d6m]]


CForm[FullSimplify[d7Z]]


CForm[FullSimplify[d7m]]


FullSimplify[d5m*kappa]
FullSimplify[d5m*kappa]
FullSimplify[d6Z*kappa]
FullSimplify[d6Z*kappa]
FullSimplify[d7Z*kappa]
FullSimplify[d7Z*kappa]







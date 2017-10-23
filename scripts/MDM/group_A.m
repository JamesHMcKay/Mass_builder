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

SEtotal=0;
SEtotalc=0;
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_23_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{24,25}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_26_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{27,28}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{29,30}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{31,34}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{32,33}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{35,36,37,38}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{39,40}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_41_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{42,43}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 


 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{11,12,15,16}_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F5_{13,14,17,18}_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 

Csum1=makeFiniteAmplitude[SEtotal,-1,D];
Csum2=makeFiniteAmplitude[SEtotal,-2,D];

Cct1 = makeFiniteAmplitude[SEtotalc,-1,D];
Cct2 = makeFiniteAmplitude[SEtotalc,-2,D];

SEtotal=0;
SEtotalc=0;

(* Doubly charged particle *)


 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_17_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{18,19}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_20_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{21,22}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{23,24}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_25_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_26_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_27_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_28_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_29_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{30,31}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 

 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{9,12}_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{10,13}_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F6_{11,14}_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 

CCsum1=makeFiniteAmplitude[SEtotal,-1,D];
CCsum2=makeFiniteAmplitude[SEtotal,-2,D];

CCct1 = makeFiniteAmplitude[SEtotalc,-1,D];
CCct2 = makeFiniteAmplitude[SEtotalc,-2,D];
(* Neutral particle *)

SEtotal=0;
SEtotalc=0;

 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{13,14,15,16}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{17,18}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{19,20}_2.mx"]; 
 SEtotal = SEtotal + SelfEnergyFinite*kappa; 

 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{5,6}_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 
 Get["/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/output/math_data_F7_{7,8}_2c.mx"]; 
 SEtotalc = SEtotalc + SelfEnergyFinite*kappa; 

Nsum1=makeFiniteAmplitude[SEtotal,-1,D];
Nsum2=makeFiniteAmplitude[SEtotal,-2,D];

Nct1 = makeFiniteAmplitude[SEtotalc,-1,D];
Nct2 = makeFiniteAmplitude[SEtotalc,-2,D];

{Nsum1,Nsum2,Csum1,Csum2,Nct1,Nct2,Cct1,Cct2,CCsum1,CCsum2,CCct1,CCct2}={Nsum1,Nsum2,Csum1,Csum2,Nct1,Nct2,Cct1,Cct2,CCsum1,CCsum2,CCct1,CCct2}/. (DiracGamma[Momentum[p, D], D].DiracGamma[6] + DiracGamma[Momentum[p, D], D].DiracGamma[7]) -> p\
/. MassBuilderP^2 -> Pair[Momentum[p],Momentum[p]]/. MassBuilderP -> Momentum[p]/. MassBuilderEpsilon -> \[Epsilon]/. MassBuilderZeta -> \[Zeta]/. MassBuilderAe[0] -> 0\
/. D -> 4 - 2 \[Epsilon]/. Pair[p, p] -> p^2/. D -> 4 - 2 \[Epsilon]/. Momentum[p] -> p\
/. MassBuilderB[a_,b_] -> TBI[4, MChi^2, {{1, a}, {1, b}}]\
/.p->MChi/. ma->0/.STW->sw/.CTW->cw;
sw=Sin[\[Theta]];
cw=Cos[\[Theta]];


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


FullSimplify[Csum2+Cct2,{sw^2+cw^2==1}]
FullSimplify[Nsum2+Nct2,{sw^2+cw^2==1}]
FullSimplify[CCsum2+CCct2/.mz->mw/cw,{sw^2+cw^2==1}]


CForm[dg1]


CForm[dg2]


CForm[dg3]


CForm[dg4]


CForm[dg5]


CForm[dg6]


CForm[dg7]


CForm[dg8]

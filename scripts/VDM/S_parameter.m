(* ::Package:: *)

Quit[]


$LoadTARCER = True;
$LoadFeynArts = True;
<< FeynCalc/FeynCalc.m;
AppendTo[$Path, 
  "/Users/jamesmckay/Documents/Programs/Mass_builder/src/"];
<< MassBuilder.m;
SetOptions[DiracSlash, Dimension -> D, 
 FeynCalcInternal -> True]; SetOptions[DiracTrace, 
 DiracTraceEvaluate -> True]; null = 0;


(* Determine one-loop counter-term couplings *)

path = "/Users/jamesmckay/Documents/Programs/Mass_builder/";
kappa = 1/(16 \[Pi]^2);
sw=Sin[\[Theta]];
cw=Cos[\[Theta]];
SW=Sin[\[Theta]];
CW=Cos[\[Theta]];

V1SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/VDM/output/math_data_V1_19_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;


Get[FileNameJoin[{path, "/models/VDM/output/math_data_V1_3_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

SE1temp = makeFiniteAmplitude[V1SE*kappa, 0, D];

V2SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/VDM/output/math_data_V2_5_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;


Get[FileNameJoin[{path, "/models/VDM/output/math_data_V2_26_1.mx"}]];
V2SE = V2SE + SelfEnergyFinite*kappa;

SE2temp = makeFiniteAmplitude[V2SE*kappa, 0, D];
V3SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/VDM/output/math_data_V1_V2_3_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/VDM/output/math_data_V1_V2_19_1.mx"}]];
V3SE = V3SE + SelfEnergyFinite*kappa;

SE3temp = makeFiniteAmplitude[V3SE*kappa, 0, D];


{SE1,SE2,SE3} = FullSimplify[{SE1temp,SE2temp,SE3temp}] /. MassBuilderEpsilon -> 0 /. Pair[Momentum[p], Momentum[p]] -> p^2;
{SE1,SE2,SE3} = {SE1,SE2,SE3} /. TBI[4, p^2, {{1, mw}, {1, MVp}}] -> I*Bwc[s]\
/. TBI[4, p^2, {{1, MVp}, {1, MVp}}] -> I*Bcc[s]\
/. TAI[4, 0, {{1, mw}}] -> I*Aw \
/. TAI[4, 0, {{1, MVp}}] -> I*Ac\
/. TAI[4, 0, {{1, ma}}] -> I*Aa \
/. TAI[4, 0, {{1, mz}}] -> I*Az;

{SE1,SE2,SE3}={SE1,SE2,SE3}/.p^2->s/.p->Sqrt[s];
{D1,D2,D3} = D[{SE1,SE2,SE3}, s];
{D1final,D2final,D3final}=Simplify[{D1,D2,D3}\
/.Derivative[1][Bzc][s]->dBzc\
/.Derivative[1][Bcc][s]->dBcc\
/.Derivative[1][Bcz][s]->dBzc\
/.Derivative[1][Bwc][s]->dBwc\
/.Derivative[1][Bcw][s]->dBwc\
/.Derivative[1][Bac][s]->dBac\
/.Derivative[1][Bca][s]->dBac\
/.Bzc[s]->Bcz/.Bca[s]->Bca\
/.Bac[s]->Bca/.Bwc[s]->Bcw\
/.Bcw[s]->Bcw/.Bcz[s]->Bcz]/.Sqrt[s]->0/.s->0;

(*CForm[Simplify[D3final,{MVp>0}]]*)
D1final
D2final
D3final
S =FullSimplify[ ((16*Pi*cw^2*sw^2)/e^2) * ( D2final - ((cw^2-sw^2)/(cw*sw) ) * D3final-D1final)]





SE1
SE2
SE3


16^2

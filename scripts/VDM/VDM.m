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
kappa=1/(16\[Pi]^2);
(*STW=Sin[\[Theta]];
CTW=Cos[\[Theta]];*)

(* V5 *)

V5SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V5_{1,2,3,4,5,6,7,8}_1.mx"}]];
V5SE = V5SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V5_{9,10,11,12}_1.mx"}]];
V5SE = V5SE + SelfEnergyFinite*kappa;

SEc = makeFiniteAmplitude[V5SE,0,D]/.MVp->MChi/.MV0->MChi;

SEc=SEc/.TBI[4,p^2,{{1,MZ},{1,MChi}}]-> I*Bzp\
/.TBI[4,p^2,{{1,MChi},{1,0}}]->I Bap\
/.TBI[4,p^2,{{1,MW},{1,MChi}}]->I Bwp\
/.TAI[4,0,{{1,MW}}]-> I*Aw\
/.TAI[4,0,{{1,MChi}}]-> I*Ap\
/.TAI[4,0,{{1,MZ}}]-> I*Az;



(* V6 *)

V6SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V6_{1,2,3,4,5}_1.mx"}]];
V6SE = V6SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MDM/output/math_data_V5_{6,7,8}_1.mx"}]];
V6SE = V6SE + SelfEnergyFinite*kappa;

SEn = makeFiniteAmplitude[V6SE,0,D]/.MVp->MChi/.MV0->MChi;


SEn=SEn/.TBI[4,p^2,{{1,MZ},{1,MChi}}]-> I*Bzp\
/.TBI[4,p^2,{{1,MChi},{1,0}}]->I Bap\
/.TBI[4,p^2,{{1,MW},{1,MChi}}]->I Bwp\
/.TAI[4,0,{{1,MW}}]-> I*Aw\
/.TAI[4,0,{{1,MChi}}]-> I*Ap\
/.TAI[4,0,{{1,MZ}}]-> I*Az;






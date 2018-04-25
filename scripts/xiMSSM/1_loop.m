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
(*sw=Sin[\[Theta]];
cw=Cos[\[Theta]];*)
(*F11*)

F11SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F11_g1_1_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F11_g1_2_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F11_g1_3_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

F11SEdiv = makeFiniteAmplitude[F11SE, -1, D];

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F11_g1_1_1c.mx"}]];

F11ct = makeFiniteAmplitude[SelfEnergyFinite*kappa, -1, D] /. 
    MassBuilderCTM1 -> 0 /. MassBuilderCTZ1 -> 0;

eq1 = FullSimplify[Coefficient[F11ct + F11SEdiv, p]];
eq2 = FullSimplify[Coefficient[F11ct + F11SEdiv, p, 0]];
solF11 = Solve[{eq1 == 0, eq2 == 0}, {d2Z, d2m}]

(*F12*)

F12SE = 0;

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F12_g1_1_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F12_g1_2_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F12_g1_3_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

F12SEdiv = makeFiniteAmplitude[F12SE, -1, D];
Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F12_g1_1_1c.mx"}]];
F12ct = makeFiniteAmplitude[SelfEnergyFinite*kappa, -1, D] /. 
    MassBuilderCTM1 -> 0 /. MassBuilderCTZ1 -> 0;
eq1 = FullSimplify[Coefficient[F12ct + F12SEdiv, p]];
eq2 = FullSimplify[Coefficient[F12ct + F12SEdiv, p, 0]];
solF12 = Solve[{eq1 == 0, eq2 == 0}, {d1Z, d1m}]
Set @@@ solF11[[1]];
Set @@@ solF12[[1]];


Clear[d1Z,d2Z,d1m,d2m]


(* Determine one-loop counter-term couplings *)

path = "/Users/jamesmckay/Documents/Programs/Mass_builder/";
kappa = 1/(16 \[Pi]^2);
(*sw=Sin[\[Theta]];
cw=Cos[\[Theta]];*)
(*F11*)

F11SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F11_g1_{1,2}_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

F11SEdiv = makeFiniteAmplitude[F11SE, -1, D];

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F11_g1_1_1c.mx"}]];

F11ct = makeFiniteAmplitude[SelfEnergyFinite*kappa, -1, D] /. 
    MassBuilderCTM1 -> 0 /. MassBuilderCTZ1 -> 0;

eq1 = FullSimplify[Coefficient[F11ct + F11SEdiv, p]];
eq2 = FullSimplify[Coefficient[F11ct + F11SEdiv, p, 0]];
solF11 = Solve[{eq1 == 0, eq2 == 0}, {d2Z, d2m}]

(*F12*)

F12SE = 0;

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F12_g1_{1,2,3}_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

F12SEdiv = makeFiniteAmplitude[F12SE, -1, D];
Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F12_g1_1_1c.mx"}]];
F12ct = makeFiniteAmplitude[SelfEnergyFinite*kappa, -1, D] /. 
    MassBuilderCTM1 -> 0 /. MassBuilderCTZ1 -> 0;
eq1 = FullSimplify[Coefficient[F12ct + F12SEdiv, p]];
eq2 = FullSimplify[Coefficient[F12ct + F12SEdiv, p, 0]];
solF12 = Solve[{eq1 == 0, eq2 == 0}, {d1Z, d1m}]
Set @@@ solF11[[1]];
Set @@@ solF12[[1]];


F11SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F11_g1_1_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F11_g1_2_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F11_g1_3_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F11_g1_1_1c.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;
F11SE1 = F11SE /. MassBuilderCTZ1 -> 0 /. MassBuilderCTZ2 -> 0 /. 
     MassBuilderCTM1 -> 0 /. MassBuilderCTM2 -> 0 /. 
   D -> 4 - 2 MassBuilderEpsilon;

F11SE1 = makeFiniteAmplitude[F11SE,0,D];

(* F12 *)
F12SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/xiMSSM///output/math_data_F12_g1_1_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/xiMSSM///output/math_data_F12_g1_2_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/xiMSSM///output/math_data_F12_g1_3_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/xiMSSM///output/math_data_F12_g1_1_1c.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

F12SE = F12SE/.MassBuilderCTZ1->0/.MassBuilderCTZ2->0/.MassBuilderCTM1->0/.MassBuilderCTM2->0/.D->4-2MassBuilderEpsilon;

F12SE1 = makeFiniteAmplitude[F12SE,0,D];


F11SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F11_g1_{1,2}_1.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, 
    "/models/xiMSSM///output/math_data_F11_g1_1_1c.mx"}]];
F11SE = F11SE + SelfEnergyFinite*kappa;
F11SE = F11SE /. MassBuilderCTZ1 -> 0 /. MassBuilderCTZ2 -> 0 /. 
     MassBuilderCTM1 -> 0 /. MassBuilderCTM2 -> 0 /. 
   D -> 4 - 2 MassBuilderEpsilon;
F11SE2 = makeFiniteAmplitude[F11SE,0,D];

(* F12 *)
F12SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/xiMSSM///output/math_data_F12_g1_{1,2,3}_1.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/xiMSSM///output/math_data_F12_g1_1_1c.mx"}]];
F12SE = F12SE + SelfEnergyFinite*kappa;

F12SE = F12SE/.MassBuilderCTZ1->0/.MassBuilderCTZ2->0/.MassBuilderCTM1->0/.MassBuilderCTM2->0/.D->4-2MassBuilderEpsilon;
F12SE2 = makeFiniteAmplitude[F12SE,0,D];


FullSimplify[F11SE1-F11SE2]
FullSimplify[F12SE1-F12SE2]





FullSimplify[F12SE1-F11SE1/.Xi->1]

(* ::Package:: *)

(* ::Section:: *)
(* Compute counter-terms for massive QED from the previously computed amplitudes*)


(* ::Subsection:: *)
(* Load required packages *)


$LoadTARCER = True;
$LoadFeynArts = True;
<< FeynCalc/FeynCalc.m;
AppendTo[$Path, "/Users/jamesmckay/Documents/Programs/Mass_builder/src/"];
<<MassBuilder.m;



(* ::Subsection:: *)
(* We get these from the models/QED2/output/ directory *)
(* Change the path below to the location of your Mass_Builder/models directory *)


path = "/Users/jamesmckay/Documents/Programs/Mass_builder/";
kappa=1/(16\[Pi]^2);
Get[FileNameJoin[{path, "/models/QED2/output/math_data_V1_1_1.mx"}]];
SE0V =SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/QED2/output/math_data_V1_1_1c.mx"}]];
ct0V =SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/QED2/output/math_data_F02_g1_1_1.mx"}]];
SE0 =SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/QED2/output/math_data_F02_g1_1_1c.mx"}]];
ct0 =SelfEnergyFinite*kappa;
Get[FileNameJoin[{path, "/models/QED2/output/math_data_F02_g1_1_2.mx"}]];
SE3 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/QED2/output/math_data_F02_g1_2_2.mx"}]];
SE2 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/QED2/output/math_data_F02_g1_5_2.mx"}]];
SE5 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/QED2/output/math_data_F02_g1_1_2c.mx"}]];
ct2 = SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/QED2/output/math_data_F02_g1_2_2c.mx"}]];
ct3 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/QED2/output/math_data_F02_g1_3_2c.mx"}]];
ct1 = 2 * SelfEnergyFinite*kappa^2;

{SE0V,ct0V,SE0,SE2,SE3,SE5,ct0,ct1,ct2,ct3} = {SE0V,ct0V,SE0,SE2,SE3,SE5,ct0,ct1,ct2,ct3}/. (DiracGamma[Momentum[p, D], D].DiracGamma[6] + DiracGamma[Momentum[p, D], D].DiracGamma[7]) -> p\
/. MassBuilderP^2 -> p^2/. MassBuilderP -> Momentum[p]/. MassBuilderEpsilon -> \[Epsilon]/. MassBuilderZeta -> \[Zeta]/. MassBuilderAe[0] -> 0/. Pair[Momentum[p], Momentum[p]] -> p^2/. Momentum[p, D] -> p\
/. D -> 4 - 2 \[Epsilon]/. Pair[p, p] -> p^2/. D -> 4 - 2 \[Epsilon]/. Momentum[p] -> p/. MassBuilderB[ma, mf] -> TBI[4, p^2, {{1, mf}, {1, ma}}] /.  MassBuilderB[mf, ma] -> TBI[4, p^2, {{1, mf}, {1, ma}}];


(* ::Subsection:: *)
(* Obtain the one and two-loop counter-terms *)


eq1 = Coefficient[ Coefficient[ct0, \[Epsilon], -1] +  Coefficient[SE0, \[Epsilon], -1], p, 0] /. MassBuilderCTM1 -> 0;
eq2 = Coefficient[ Coefficient[ct0, \[Epsilon], -1] +  Coefficient[SE0, \[Epsilon], -1], p, 1] /. MassBuilderCTZ1 -> 0;
Print["We take Z1 = -EL^2"]
Print["The other counter-term couplings are calculated to be:"]
sol1 = Solve[{eq1 == 0, eq2 == 0}, {Z2z, Z2m}]
Set @@@ sol1[[1]];
eq1V = Coefficient[Coefficient[Series[ct0V, {\[Epsilon], 0, 0}], \[Epsilon], -1] + Coefficient[Series[SE0V, {\[Epsilon], 0, 0}], \[Epsilon], -1], p,0] /. MassBuilderCTM1 -> 0;
eq2V = Coefficient[Coefficient[Series[ct0V, {\[Epsilon], 0, 0}], \[Epsilon], -1] + Coefficient[Series[SE0V, {\[Epsilon], 0, 0}], \[Epsilon], -1], p,2] /. MassBuilderCTZ1 -> 0;
sol1 = Solve[{eq1V == 0, eq2V == 0}, {Z3, Z3m}]
Set @@@ sol1[[1]];
div1 = FullSimplify[Series[ct1, {\[Epsilon], 0, 0}] + Series[SE3, {\[Epsilon], 0, -1}] /. Z1 -> -EL^2] ;
div2 = FullSimplify[ Series[ct2, {\[Epsilon], 0, -1}] + Series[SE5, {\[Epsilon], 0, -1}]] ;
div3 = FullSimplify[Series[ct3, {\[Epsilon], 0, -1}] +  Series[SE2, {\[Epsilon], 0, -1}]];
div4 = FullSimplify[ Series[ct0, {\[Epsilon], 0, -1}] + Series[SE0, {\[Epsilon], 0, -1}]];
eq21 = FullSimplify[Coefficient[ Coefficient[div1 + div2 + div3 + div4, p, 1]/kappa^2, \[Epsilon], -1]];
eq22 = FullSimplify[ Coefficient[Coefficient[div1 + div2 + div3 + div4, p, 1]/kappa^2, \[Epsilon], -2]];
eq11 = FullSimplify[ Coefficient[ Coefficient[div1 + div2 + div3 + div4, p, 0]/kappa^2 /. (DiracGamma[p, 4].DiracGamma[6] + DiracGamma[p, 4].DiracGamma[7]) -> p, \[Epsilon], -1]];
eq12 = FullSimplify[ Coefficient[Coefficient[div1 + div2 + div3 + div4, p, 0]/ kappa^2 /. (DiracGamma[p, 4].DiracGamma[6] + DiracGamma[p, 4].DiracGamma[7]) -> p, \[Epsilon], -2]];
sol2 = FullSimplify[Solve[{eq11 == 0, eq21 == 0}, {MassBuilderCTZ1, MassBuilderCTM1}]]
sol3 = FullSimplify[Solve[{eq12 == 0, eq22 == 0}, {MassBuilderCTZ2, MassBuilderCTM2}]]
Set @@@ sol2[[1]];
Set @@@ sol3[[1]];



(* ::Subsection:: *)
(* Check that the total amplitude is now divergence free *)


Print["Checking that result is divergence free, coefficients of \
1/epsilon and 1/epsilon^2 are:"]
SE = SE0 + SE2 + SE3 + SE5 + ct1 + ct3 + ct2 + ct0/. (DiracGamma[Momentum[p, D], D].DiracGamma[6] + DiracGamma[Momentum[p, D], D].DiracGamma[7]) -> p\
/. MassBuilderP^2 -> p^2/. MassBuilderP -> Momentum[p]/. MassBuilderEpsilon -> \[Epsilon]/. MassBuilderZeta -> \[Zeta]/. MassBuilderAe[0] -> 0/. Pair[Momentum[p], Momentum[p]] -> p^2\
/. Momentum[p, D] -> p/. D -> 4 - 2 \[Epsilon]/. Pair[p, p] -> p^2/. D -> 4 - 2 \[Epsilon]/. Momentum[p] -> p/.MassBuilderB[ma, mf] -> TBI[4, p^2, {{1, mf}, {1, ma}}]/.MassBuilderB[mf, ma] -> TBI[4, p^2, {{1, mf}, {1, ma}}];
SE = Series[SE, {\[Epsilon], 0, -1}] /. Z1 -> -EL^2;
eq1 = FullSimplify[Coefficient[Series[SE, {\[Epsilon], 0, 0}], \[Epsilon], -1]]
eq1 = FullSimplify[Coefficient[Series[SE, {\[Epsilon], 0, 0}], \[Epsilon], -2]]




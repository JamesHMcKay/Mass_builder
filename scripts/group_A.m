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

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_17_2.mx"}]];
SE13 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_18_2.mx"}]];
SE14 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_19_2.mx"}]];
SE15 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_20_2.mx"}]];
SE16 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_21_2.mx"}]];
SE17 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_22_2.mx"}]];
SE18 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_23_2.mx"}]];
SE19 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_24_2.mx"}]];
SE20 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_25_2.mx"}]];
SE21 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_26_2.mx"}]];
SE22 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_27_2.mx"}]];
SE23 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_28_2.mx"}]];
SE24 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_29_2.mx"}]];
SE25 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_30_2.mx"}]];
SE26 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_31_2.mx"}]];
SE27 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_7_2c.mx"}]];
SE28 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_8_2c.mx"}]];
SE29 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_9_2c.mx"}]];
SE30 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_10_2c.mx"}]];
SE31 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_11_2c.mx"}]];
SE32 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F11_g1_12_2c.mx"}]];
SE33 =SelfEnergyFinite*kappa^2;



(* Particle 2 *)


Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_17_2.mx"}]];
SE1 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_18_2.mx"}]];
SE2 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_19_2.mx"}]];
SE3 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_20_2.mx"}]];
SE4 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_21_2.mx"}]];
SE5 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_22_2.mx"}]];
SE6 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_23_2.mx"}]];
SE7 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_24_2.mx"}]];
SE8 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_25_2.mx"}]];
SE9 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_26_2.mx"}]];
SE10 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_27_2.mx"}]];
SE11 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_28_2.mx"}]];
SE12 =SelfEnergyFinite*kappa^2;




Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_9_2c.mx"}]];
SE34 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_10_2c.mx"}]];
SE35 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_11_2c.mx"}]];
SE36 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_12_2c.mx"}]];
SE37 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_13_2c.mx"}]];
SE38 =SelfEnergyFinite*kappa^2;

Get[FileNameJoin[{path, "/models/massiveMSSM/output/math_data_F12_g1_14_2c.mx"}]];
SE39 =SelfEnergyFinite*kappa^2;





sum1=makeFiniteAmplitude[SE1+SE2+SE3+SE4+SE5+SE6+SE7+SE8+SE9+SE10+SE11+SE12+SE34+SE35+SE36+SE37+SE38+SE39,-2,D];
sum2=makeFiniteAmplitude[SE13+SE14+SE15+SE16+SE17+SE18+SE19+SE20+SE21+SE22+SE23+SE24+SE25+SE26+SE27+SE28+SE29+SE30+SE31+SE32+SE33,-2,D];


{sum1,sum2}={sum1,sum2}/. (DiracGamma[Momentum[p, D], D].DiracGamma[6] + DiracGamma[Momentum[p, D], D].DiracGamma[7]) -> p\
/. MassBuilderP^2 -> Pair[Momentum[p],Momentum[p]]/. MassBuilderP -> Momentum[p]/. MassBuilderEpsilon -> \[Epsilon]/. MassBuilderZeta -> \[Zeta]/. MassBuilderAe[0] -> 0\
/. D -> 4 - 2 \[Epsilon]/. Pair[p, p] -> p^2/. D -> 4 - 2 \[Epsilon]/. Momentum[p] -> p\
/. MassBuilderB[MChi, mz] -> TBI[4, MChi^2, {{1, MChi}, {1, mz}}]\
/. MassBuilderB[mz, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, mz}}]\
/. MassBuilderB[MChi, mw] -> TBI[4, MChi^2, {{1, MChi}, {1, mw}}]\
/. MassBuilderB[mw, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, mw}}]\
/. MassBuilderB[0, MChi] -> TBI[4, MChi^2, {{1, MChi}, {1, 0}}]\
/. MassBuilderB[MChi, 0] -> TBI[4, MChi^2, {{1, MChi}, {1, 0}}]\
/.p->MChi/. ma->0/.STW->sw/.CTW->cw;


eq2=Coefficient[FullSimplify[sum1,{sw^2+cw^2==1,Cw1+Cw2==-4}],TBI[4, MChi^2, {{1, MChi}, {1, 0}}]]
eq3=Coefficient[FullSimplify[sum1,{sw^2+cw^2==1,Cw1+Cw2==-4}],TBI[4, MChi^2, {{1, MChi}, {1, mw}}]]
eq4=Coefficient[FullSimplify[sum1,{sw^2+cw^2==1,Cw1+Cw2==-4}],TBI[4, MChi^2, {{1, MChi}, {1, mz}}]]
eq5=Coefficient[FullSimplify[sum1,{sw^2+cw^2==1,Cw1+Cw2==-4}],TAI[4, 0, {1, MChi}]]
eq6=Coefficient[FullSimplify[sum1,{sw^2+cw^2==1,Cw1+Cw2==-4}],TAI[4, 0, {1, mw}]]
eq7=Coefficient[FullSimplify[sum1,{sw^2+cw^2==1,Cw1+Cw2==-4}],TAI[4, 0, {1, mz}]]


Ca=-2 sw;
Cw1+Cw2==-4;
Cz2=-2cw;
sw=Sin[Theta];
cw=Cos[Theta];
Solve[eq2==0,{dg2}]
Solve[eq3==0,{dg2}]
Solve[eq4==0,{dg2}]
Solve[eq5==0,{dg2}]
Solve[eq6==0,{dg2}]
Solve[eq7==0,{dg2}]


eq1=FullSimplify[sum1-sum2/.dg2->4*e^3/sw^3,{sw^2+cw^2==1,Cw1+Cw2==-4}]

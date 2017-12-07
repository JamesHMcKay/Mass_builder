(* ::Package:: *)

Quit[]


(* ::Section:: *)
(**)


$LoadAddOns = {"FeynHelpers"}; 
 $LoadTARCER = True; 
 $LoadFeynArts = True; 
 << FeynCalc/FeynCalc.m; 
 AppendTo[$Path, "/Users/jamesmckay/Documents/Programs/Mass_builder/src/"]; 
 << MassBuilder.m; 
 
   


(* ::Section:: *)
(**)


SetOptions[DiracSlash, Dimension -> D, FeynCalcInternal -> True]; 
 SetOptions[DiracTrace, DiracTraceEvaluate -> True]; 
 null=0; 
 MassBuilderA[mass_, D_] := TAI[D, 0, {{1, mass}}]; 
 MassBuilderB[mass1_, mass2_, D_] := TBI[D, Pair[Momentum[p, D],Momentum[p, D]], {{1, mass1}, {1, mass2}}]; 
 MassBuilderJ[mass1_, mass2_, mass3_, D_] := TJI[D, Pair[Momentum[p, D],Momentum[p, D]], {{1, mass1}, {1, mass2}, {1, mass3}}]; 
 MassBuilderT[mass1_, mass2_, mass3_, D_] := TJI[D, Pair[Momentum[p, D],Momentum[p, D]], {{2, mass1}, {1, mass2}, {1, mass3}}]; 
 MassBuilderK[mass1_, mass2_, mass3_, D_] := TJI[D, 0, {{1, mass1}, {1, mass2}, {1, mass3}}]; 
 MassBuilderV[mass1_, mass2_, mass3_, mass4_, D_] := TVI[D, Pair[Momentum[p, D],Momentum[p, D]], {{1, mass1}, {1, mass2}, {1, mass3}, {1, mass4}}]; 
 MassBuilderF[mass1_, mass2_, mass3_, mass4_, mass5_, D_] := TFI[D, Pair[Momentum[p, D],Momentum[p, D]], {{1, mass1}, {1, mass2}, {1, mass3}, {1,mass4}, {1,mass5}}]; 
 FCGV["EL"] = e; 
 FCGV["SW"] = sw; 
 FCGV["CW"] = cw; 
 FCGV["MW"] = mw; 
 FCGV["MZ"] = mz; 
 FCGV["ME"] = me; 
 FCGV["MM"] = mm; 
 FCGV["ML"] = ml; 
 FCGV["MU"] = mu; 
 FCGV["MT"] = mt; 
 FCGV["MD"] = md; 
 FCGV["MS"] = ms; 
 FCGV["MB"] = mb; 
 FCGV["MF"] = mf; 
 ZNeu[1,1] = 0; 
 ZNeu[1,2] = 1; 
 ZNeu[1,3] = 0; 
 ZNeu[1,4] = 0; 
 ZNeu[2,1] = 0; 
 ZNeu[2,2] = 0; 
 ZNeu[2,3] = 0; 
 ZNeu[2,4] = 0; 
 ZNeu[3,1] = 0; 
 ZNeu[3,2] = 0; 
 ZNeu[3,3] = 0; 
 ZNeu[3,4] = 0; 
 ZNeu[4,1] = 0; 
 ZNeu[4,2] = 0; 
 ZNeu[4,3] = 0; 
 ZNeu[4,4] = 0; 
 UCha[1,1] = 1; 
 UCha[1,2] = 0; 
 UCha[2,1] = 0; 
 UCha[2,2] = 0; 
 VCha[1,1] = 1; 
 VCha[1,2] = 0; 
 VCha[2,1] = 0; 
 VCha[2,2] = 0; 
 CKM = IndexDelta; 
 MCha[1] = MChi; 
 MCha[2] = MChi; 
 MNeu[1] = MChi; 
 MNeu[2] = MChi; 
 Mh0 = mh; 
 SBA = 1; 
 t12 = CreateTopologies[1, 1 -> 1, ExcludeTopologies -> Internal]; 
 alldiags = InsertFields[t12, {F[11,{1}]} -> {F[11,{1}]},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Restrictions -> {WinoLimit,WinoCouplings},Model -> "/Users/jamesmckay/Documents/Programs/Mass_builder/models/xiMSSM/xiMSSM"]; 
 subdiags0 =   DiagramExtract[alldiags, {1,2}]; 
 amp0 = FCFAConvert[CreateFeynAmp[subdiags0,Truncated -> True, GaugeRules -> {}], IncomingMomenta -> {p}, OutgoingMomenta -> {p}, LoopMomenta -> {k1, k2} ,UndoChiralSplittings -> True,TransversePolarizationVectors -> {p}, DropSumOver -> True, List -> False,ChangeDimension -> D] // Contract // FCTraceFactor; 
 amp0 = amp0 /. GaugeXi[Z] -> Xi^2 /. GaugeXi[A] -> Xi^2 /. GaugeXi[W] -> Xi^2 /. GaugeXi[P] -> Xi^2 /. GaugeXi[Wp] -> Xi^2; 
  amp0 = amp0*16*Pi^4; 
 ToString[amp0/.Sqrt[Xi^2]->Xi,InputForm]
masses = List[MChi,mw,mw*Xi]; 
 massesExpand = List[1,1,1]; 
 
amp0 = amp0 /. Index[Generation, 1] -> 1; 
 amp0 = amp0 /. Index[Generation, 2] -> 2; 
 amp0 = amp0 /. Index[Generation, 3] -> 3; 
  ampsSE1 = (amp0) // DiracSimplify // TID[#, k1] & // DiracSimplify; 
 amp1 = ampsSE1 //FDS[#,k1,k2]&; 
 SetOptions[Eps, Dimension -> D]; 
  fullamp0 = (amp1) // DiracSimplify // TID[#, k1] & // DiracSimplify; 
 tfiamp0 = fullamp0 // ToTFI[#, k1, p] & // ChangeDimension[#, D] &; 
 
SE = TarcerRecurse[tfiamp0]; 
 SEk = (1/(4 Pair[Momentum[p, D],Momentum[p, D]])) DiracTrace[ DiracGamma[Momentum[p, D], D].SE ]; 
 SEm = (1/4) DiracTrace[ SE ]; 
 SE = p*SEk+SEm; 
 SE = SE /. Pair[Momentum[Polarization[p, -I, Transversality -> True], D], Momentum[Polarization[p, I, Transversality -> True], D]] -> -1 ; 
 SE = SE /. Sqrt[Xi^2]->Xi; 
 
SelfEnergyFinite = expandBasisIntegrals[SE, masses, massesExpand, MassBuilderA,MassBuilderB, MassBuilderJ, MassBuilderK, MassBuilderT, MassBuilderV, MassBuilderF]; 
 DumpSave["/Users/jamesmckay/Documents/Programs/Mass_builder/models/xiMSSM/output/math_data_F11_g1_{1,2}_1.mx", SelfEnergyFinite]; 
 
SelfEnergyFinite = makeFiniteAmplitude[SelfEnergyFinite,0, D]; 
 SelfEnergyFinite = SelfEnergyFinite/.(Pair[LorentzIndex[Lor1], Momentum[p]]*Pair[LorentzIndex[Lor2], Momentum[p]] -Pair[LorentzIndex[Lor1], LorentzIndex[Lor2]]*Pair[Momentum[p], Momentum[p]]) -> Pair[Momentum[p],Momentum[p]] ; 
 SelfEnergyFinite = SelfEnergyFinite/.(-Pair[LorentzIndex[Lor1], Momentum[p]]*Pair[LorentzIndex[Lor2], Momentum[p]] +Pair[LorentzIndex[Lor1], LorentzIndex[Lor2]]*Pair[Momentum[p], Momentum[p]]) -> -Pair[Momentum[p],Momentum[p]] ; 
 SelfEnergyFinite = Simplify[SelfEnergyFinite/.MassBuilderP^2 -> Pair[Momentum[p],Momentum[p]] /. MassBuilderP -> p /. MassBuilderQ2->Q2 /. MassBuilderZeta-> Zeta,TimeConstraint->100000]; 
 massesSmall = List[0,0,0]; 
 SelfEnergyFinite = Simplify[implementTbar[SelfEnergyFinite,masses, massesSmall,MassBuilderA, MassBuilderB, MassBuilderT],TimeConstraint->100000]; 
 



 << /Users/jamesmckay/Documents/Programs/Mass_builder/output/math_1_0.m; 
 
<< /Users/jamesmckay/Documents/Programs/Mass_builder/output/math_2_0.m; 
 
<< /Users/jamesmckay/Documents/Programs/Mass_builder/output/math_3_0.m; 
 


sw=Sin[\[Theta]];
cw=Cos[\[Theta]];
Xi = \[Xi];


Coefficient[CAc*I/.e->g*sw,p,1]
Coefficient[CAc*I,p,0]


Coefficient[CAw*I/.e->g*sw ,p,1]
Coefficient[CAw*I ,p,0]


Coefficient[CAwx*I/.e->g*sw,p,1]
Coefficient[CAwx*I,p,0]


Coefficient[CBcw*I/.e->g*sw ,p,1]
Coefficient[CBcw*I/.e->g*sw ,p,0]


Coefficient[CBcwx*I/.e->g*sw ,p,1]
Coefficient[CBcwx*I/.e->g*sw ,p,0]


Coefficient[diff/.e->g*sw,p,1]
Coefficient[diff/.e->g*sw,p,0]




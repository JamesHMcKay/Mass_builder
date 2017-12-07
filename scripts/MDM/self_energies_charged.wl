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
 ZUL[1,1] = 1; 
 ZUL[1,2] = 0; 
 ZUL[1,3] = 0; 
 ZUL[2,1] = 0; 
 ZUL[2,2] = 1; 
 ZUL[2,3] = 0; 
 ZUL[3,1] = 0; 
 ZUL[3,2] = 0; 
 ZUL[3,3] = 1; 
 ZDL[1,1] = 1; 
 ZDL[1,2] = 0; 
 ZDL[1,3] = 0; 
 ZDL[2,1] = 0; 
 ZDL[2,2] = 1; 
 ZDL[2,3] = 0; 
 ZDL[3,1] = 0; 
 ZDL[3,2] = 0; 
 ZDL[3,3] = 1; 
 ZEL[1,1] = 1; 
 ZEL[1,2] = 0; 
 ZEL[1,3] = 0; 
 ZEL[2,1] = 0; 
 ZEL[2,2] = 1; 
 ZEL[2,3] = 0; 
 ZEL[3,1] = 0; 
 ZEL[3,2] = 0; 
 ZEL[3,3] = 1; 
 t12 = CreateTopologies[1, 1 -> 1, ExcludeTopologies -> Internal]; 
 alldiags = InsertFields[t12, {F[5]} -> {F[5]},InsertionLevel -> {Particles}, GenericModel -> Lorentz,Restrictions -> {},Model -> "/Users/jamesmckay/Documents/Programs/Mass_builder/models/MDM/MDM"]; 
 subdiags0 =   DiagramExtract[alldiags, {1,2,3,4}]; 
 amp0 = FCFAConvert[CreateFeynAmp[subdiags0,Truncated -> True], IncomingMomenta -> {p}, OutgoingMomenta -> {p}, LoopMomenta -> {k1, k2} ,UndoChiralSplittings -> True,TransversePolarizationVectors -> {p}, DropSumOver -> False, List -> False,ChangeDimension -> D] // Contract // FCTraceFactor; 
 amp0 = amp0 /.SumOver[SUNFIndex[x_], 3]->3; 
 amp0 = amp0*16*Pi^4; 
 ToString[amp0/.Sqrt[Xi^2]->Xi,InputForm]
masses = List[MChi,ma,mw,mz]; 
 massesExpand = List[1,1,1,1]; 
 
amp0 = amp0 /. Index[Generation, 1] -> 1; 
 amp0 = amp0 /. Index[Generation, 2] -> 2; 
 amp0 = amp0 /. Index[Generation, 3] -> 3; 
  ampsSE1 = (amp0) // DiracSimplify // TID[#, k1] & // DiracSimplify; 
 ampsSE3 = FIREBurn[ampsSE1, {k1}, {p}] // FDS[#, k1] &; 
 ampsSE4 = ampsSE3 // Collect2[#, {FeynAmpDenominator}, Factoring -> Simplify] &; 
 resSE = Cancel[ampsSE4]; 
 tfiamp0 = resSE // ToTFI[#, k1, p] &; 
 
SE = TarcerRecurse[tfiamp0]; 
 SEk = (1/(4 Pair[Momentum[p, D],Momentum[p, D]])) DiracTrace[ DiracGamma[Momentum[p, D], D].SE ]; 
 SEm = (1/4) DiracTrace[ SE ]; 
 SE = p*SEk+SEm; 
 SE = SE /. Pair[Momentum[Polarization[p, -I, Transversality -> True], D], Momentum[Polarization[p, I, Transversality -> True], D]] -> -1 ; 
 SE = SE /. Sqrt[Xi^2]->Xi; 
 
SelfEnergyFinite = expandBasisIntegrals[SE, masses, massesExpand, MassBuilderA,MassBuilderB, MassBuilderJ, MassBuilderK, MassBuilderT, MassBuilderV, MassBuilderF]; 
 
SelfEnergyFinite = makeFiniteAmplitude[SelfEnergyFinite,0, D]; 
 SelfEnergyFinite = SelfEnergyFinite/.(Pair[LorentzIndex[Lor1], Momentum[p]]*Pair[LorentzIndex[Lor2], Momentum[p]] -Pair[LorentzIndex[Lor1], LorentzIndex[Lor2]]*Pair[Momentum[p], Momentum[p]]) -> Pair[Momentum[p],Momentum[p]] ; 
 SelfEnergyFinite = SelfEnergyFinite/.(-Pair[LorentzIndex[Lor1], Momentum[p]]*Pair[LorentzIndex[Lor2], Momentum[p]] +Pair[LorentzIndex[Lor1], LorentzIndex[Lor2]]*Pair[Momentum[p], Momentum[p]]) -> -Pair[Momentum[p],Momentum[p]] ; 
 SelfEnergyFinite = Simplify[SelfEnergyFinite/.MassBuilderP^2 -> Pair[Momentum[p],Momentum[p]] /. MassBuilderP -> p /. MassBuilderQ2->Q2 /. MassBuilderZeta-> Zeta,TimeConstraint->100000]; 
 massesSmall = List[0,1,0,0]; 
 SelfEnergyFinite = Simplify[implementTbar[SelfEnergyFinite,masses, massesSmall,MassBuilderA, MassBuilderB, MassBuilderT],TimeConstraint->100000]; 
 
amp1 = Simplify[Coefficient[SelfEnergyFinite/. MassBuilderA[ma,4]-> (ma^2 MassBuilderA[ma,4]),ma,-2],TimeConstraint->100000]; 
 amp2 = Simplify[SelfEnergyFinite - amp1/ma^2,TimeConstraint->100000]; 
 


 << /Users/jamesmckay/Documents/Programs/Mass_builder/output/math_1_0.m; 
 
<< /Users/jamesmckay/Documents/Programs/Mass_builder/output/math_2_0.m; 
 
<< /Users/jamesmckay/Documents/Programs/Mass_builder/output/math_3_0.m; 
 


STW=Sin[\[Theta]];
CTW=Cos[\[Theta]];


Coefficient[FullSimplify[C2AC]*I/.e->g*sw,p,1]
Coefficient[FullSimplify[C2AC]*I,p,0]


Coefficient[FullSimplify[C2Aw]*I/.e->g*sw ,p,1]
Coefficient[FullSimplify[C2Aw]*I ,p,0]


Coefficient[FullSimplify[C2Aa]*I/.e->g*sw ,p,1]
Coefficient[FullSimplify[C2Aa]*I ,p,0]


Coefficient[FullSimplify[C2Az]*I/.e->g*sw ,p,1]
Coefficient[FullSimplify[C2Az]*I ,p,0]


Coefficient[FullSimplify[C2BCw]*I/.e->g*sw ,p,1]
Coefficient[FullSimplify[C2BCw]*I/.e->g*sw ,p,0]


Coefficient[FullSimplify[C2BCa]*I/.e->g*sw ,p,1]
Coefficient[FullSimplify[C2BCa]*I/.e->g*sw ,p,0]


Coefficient[FullSimplify[C2BCz]*I/.e->g*sw ,p,1]
Coefficient[FullSimplify[C2BCz]*I/.e->g*sw ,p,0]


Coefficient[FullSimplify[difference2]/.e->g*sw,p,1]
Coefficient[FullSimplify[difference2]/.e->g*sw,p,0]




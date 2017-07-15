(* ::Package:: *)

(* ::Title::  *)
(* Introduction to Mass Builder *)


(* ::Summary::  *)
(* This tutorial is designed to introduce the procedure used to obtain an amplitude *)


(* ::Section:: *)
(* Getting the amplitude from FeynArts + FeynCalc*)


(* ::Subsection:: *)
(* Load all the required packages*)
(* Change the path below to the location of your Mass_Builder/src directory *)


$LoadTARCER = True;
$LoadFeynArts = True;
<< FeynCalc/FeynCalc.m;
AppendTo[$Path, "/Users/jamesmckay/Documents/Programs/Mass_builder/src/"];
<<MassBuilder.m;



(* ::Subsection:: *)
(* Set some options and define the MassBuilder basis functions that we will need later *)


SetOptions[DiracSlash, Dimension -> D, FeynCalcInternal -> True]; SetOptions[DiracTrace,DiracTraceEvaluate -> True]; null = 0;
MassBuilderA[mass_, D_] := TAI[D, 0, {{1, mass}}]; 
MassBuilderB[mass1_, mass2_, D_] := TBI[D, Pair[Momentum[p, D],Momentum[p, D]], {{1, mass1}, {1, mass2}}];
MassBuilderJ[mass1_, mass2_, mass3_, D_] := TJI[D, Pair[Momentum[p, D],Momentum[p, D]], {{1, mass1}, {1, mass2}, {1, mass3}}];
MassBuilderT[mass1_, mass2_, mass3_, D_] := TJI[D, Pair[Momentum[p, D],Momentum[p, D]], {{2, mass1}, {1, mass2}, {1, mass3}}];
MassBuilderK[mass1_, mass2_, mass3_, D_] := TJI[D, 0, {{1, mass1}, {1, mass2}, {1, mass3}}];
MassBuilderV[mass1_, mass2_, mass3_, mass4_, D_] := TVI[D, Pair[Momentum[p, D],Momentum[p, D]], {{1, mass1}, {1, mass2}, {1, mass3}, {1, mass4}}];
MassBuilderF[mass1_, mass2_, mass3_, mass4_, D_] :=TJI[D, Pair[Momentum[p, D],Momentum[p, D]], {{1, mass1}, {1, mass2}, {1, mass3}, {1, mass4}}];




(* ::Subsection:: *)
(* Get the amplitude *)


t12 = CreateTopologies[2, 1 -> 1, ExcludeTopologies -> Internal]; 
alldiags = InsertFields[t12, {S[1]} -> {S[1]}, InsertionLevel -> {Particles},GenericModel -> Lorentz, Restrictions -> {},Model -> "/Users/jamesmckay/Documents/Programs/Mass_builder/models/Scalar/Scalar"];
subdiags0 = DiagramExtract[alldiags, 1];
Paint[%,Numbering -> None, ColumnsXRows -> 1];
amp0 =FCFAConvert[CreateFeynAmp[subdiags0, Truncated -> True],IncomingMomenta -> {p}, OutgoingMomenta -> {p},
LoopMomenta -> {k1, k2}, UndoChiralSplittings -> True, TransversePolarizationVectors -> {p}, DropSumOver -> True, List -> False, ChangeDimension -> D] // Contract // FCTraceFactor;
amp0 = amp0*256*Pi^8;
masses = List[Ms]; 
massesExpand = List[1];
amp0 = amp0 /. Index[Generation, 1] -> 1;
amp0 = amp0 /. Index[Generation, 2] -> 2;
amp0 = amp0 /. Index[Generation, 3] -> 3;
amp1 = amp0 // FDS[#, l1, l2] &;
SetOptions[Eps,Dimension -> D];
fullamp0 = (amp1) // DiracSimplify //FCMultiLoopTID[#, {k1, k2}] & // DiracSimplify;
tfiamp0 = fullamp0 // ToTFI[#, k1, k2, p] & // ChangeDimension[#, D] &;
SE = Simplify[TarcerRecurse[tfiamp0]];
SEk = (1/(4 Pair[Momentum[p, D],Momentum[p, D]])) DiracTrace[DiracGamma[Momentum[p, D], D].SE];
SEm = (1/4) DiracTrace[SE];
SE =p*SEk + SEm;
SE = SE /. Pair[Momentum[Polarization[p, -I, Transversality -> True], D],Momentum[Polarization[p, I, Transversality -> True], D]] -> -1;
SelfEnergyFinite = expandBasisIntegrals[SE, masses, massesExpand, MassBuilderA,MassBuilderB, MassBuilderJ, MassBuilderK, MassBuilderT,MassBuilderV, MassBuilderF];

SelfEnergyFinite=\[Epsilon]*makeFiniteAmplitude[SelfEnergyFinite,1,D]+makeFiniteAmplitude[SelfEnergyFinite,0,D]+(1/\[Epsilon])*makeFiniteAmplitude[SelfEnergyFinite,-1,D]+(1/\[Epsilon]^2)*makeFiniteAmplitude[SelfEnergyFinite,-2,D];
SelfEnergyFinite=FullSimplify[SelfEnergyFinite/.MassBuilderP^2->Pair[Momentum[p],Momentum[p]]/.MassBuilderP->p/.MassBuilderQ2->Q2/.MassBuilderZeta->Zeta];
SelfEnergyFinite = FullSimplify[SelfEnergyFinite /. MassBuilderEpsilon -> \[Epsilon] /. MassBuilderP -> p];
SelfEnergySeries = FullSimplify[Series[SelfEnergyFinite,{\[Epsilon],0,0}]]




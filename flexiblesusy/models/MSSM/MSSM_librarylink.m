libMSSM = FileNameJoin[{Directory[], "models", "MSSM", "MSSM_librarylink.so"}];

FSMSSMGetSettings = LibraryFunctionLoad[libMSSM, "FSMSSMGetSettings", LinkObject, LinkObject];
FSMSSMGetSMInputParameters = LibraryFunctionLoad[libMSSM, "FSMSSMGetSMInputParameters", LinkObject, LinkObject];
FSMSSMGetInputParameters = LibraryFunctionLoad[libMSSM, "FSMSSMGetInputParameters", LinkObject, LinkObject];

FSMSSMOpenHandleLib = LibraryFunctionLoad[libMSSM, "FSMSSMOpenHandle", {{Real,1}}, Integer];
FSMSSMCloseHandle = LibraryFunctionLoad[libMSSM, "FSMSSMCloseHandle", {Integer}, Void];

FSMSSMSetLib = LibraryFunctionLoad[libMSSM, "FSMSSMSet", {Integer, {Real,1}}, Void];

FSMSSMCalculateSpectrum = LibraryFunctionLoad[libMSSM, "FSMSSMCalculateSpectrum", LinkObject, LinkObject];
FSMSSMCalculateObservables = LibraryFunctionLoad[libMSSM, "FSMSSMCalculateObservables", LinkObject, LinkObject];

FSMSSMCalculateSpectrum::error = "`1`";
FSMSSMCalculateSpectrum::warning = "`1`";

FSMSSMCalculateObservables::error = "`1`";
FSMSSMCalculateObservables::warning = "`1`";

FSMSSM::info = "`1`";
FSMSSM::nonum = "Error: `1` is not a numeric input value!";
FSMSSMMessage[s_] := Message[FSMSSM::info, s];

FSMSSMCheckIsNumeric[a_?NumericQ] := a;
FSMSSMCheckIsNumeric[a_] := (Message[FSMSSM::nonum, a]; Abort[]);

fsDefaultSettings = {
      precisionGoal -> 1.*^-4,           (* FlexibleSUSY[0] *)
      maxIterations -> 0,                (* FlexibleSUSY[1] *)
      calculateStandardModelMasses -> 0, (* FlexibleSUSY[3] *)
      poleMassLoopOrder -> 2,            (* FlexibleSUSY[4] *)
      ewsbLoopOrder -> 2,                (* FlexibleSUSY[5] *)
      betaFunctionLoopOrder -> 3,        (* FlexibleSUSY[6] *)
      thresholdCorrectionsLoopOrder -> 2,(* FlexibleSUSY[7] *)
      higgs2loopCorrectionAtAs -> 1,     (* FlexibleSUSY[8] *)
      higgs2loopCorrectionAbAs -> 1,     (* FlexibleSUSY[9] *)
      higgs2loopCorrectionAtAt -> 1,     (* FlexibleSUSY[10] *)
      higgs2loopCorrectionAtauAtau -> 1, (* FlexibleSUSY[11] *)
      forceOutput -> 0,                  (* FlexibleSUSY[12] *)
      topPoleQCDCorrections -> 1,        (* FlexibleSUSY[13] *)
      betaZeroThreshold -> 1.*^-11,      (* FlexibleSUSY[14] *)
      forcePositiveMasses -> 0,          (* FlexibleSUSY[16] *)
      poleMassScale -> 0,                (* FlexibleSUSY[17] *)
      eftPoleMassScale -> 0,             (* FlexibleSUSY[18] *)
      eftMatchingScale -> 0,             (* FlexibleSUSY[19] *)
      eftMatchingLoopOrderUp -> 2,       (* FlexibleSUSY[20] *)
      eftMatchingLoopOrderDown -> 1,     (* FlexibleSUSY[21] *)
      eftHiggsIndex -> 0,                (* FlexibleSUSY[22] *)
      calculateBSMMasses -> 1,           (* FlexibleSUSY[23] *)
      parameterOutputScale -> 0          (* MODSEL[12] *)
};

fsDefaultSMParameters = {
    alphaEmMZ -> 1/127.916, (* SMINPUTS[1] *)
    GF -> 1.16637*^-5,      (* SMINPUTS[2] *)
    alphaSMZ -> 0.1184,     (* SMINPUTS[3] *)
    MZ -> 91.1876,          (* SMINPUTS[4] *)
    mbmb -> 4.18,           (* SMINPUTS[5] *)
    Mt -> 173.34,           (* SMINPUTS[6] *)
    Mtau -> 1.777,          (* SMINPUTS[7] *)
    Mv3 -> 0,               (* SMINPUTS[8] *)
    MW -> 80.385,           (* SMINPUTS[9] *)
    Me -> 0.000510998902,   (* SMINPUTS[11] *)
    Mv1 -> 0,               (* SMINPUTS[12] *)
    Mm -> 0.1056583715,     (* SMINPUTS[13] *)
    Mv2 -> 0,               (* SMINPUTS[14] *)
    md2GeV -> 0.00475,      (* SMINPUTS[21] *)
    mu2GeV -> 0.0024,       (* SMINPUTS[22] *)
    ms2GeV -> 0.104,        (* SMINPUTS[23] *)
    mcmc -> 1.27,           (* SMINPUTS[24] *)
    CKMTheta12 -> 0,
    CKMTheta13 -> 0,
    CKMTheta23 -> 0,
    CKMDelta -> 0,
    PMNSTheta12 -> 0,
    PMNSTheta13 -> 0,
    PMNSTheta23 -> 0,
    PMNSDelta -> 0,
    PMNSAlpha1 -> 0,
    PMNSAlpha2 -> 0,
    alphaEm0 -> 1/137.035999074,
    Mh -> 125.09
};

fsMSSMDefaultInputParameters = {
   TanBeta -> 0,
   SignMu -> 0,
   MassWBIN -> 0,
   Qin -> 0,
   QEWSB -> 0,
   mHd2IN -> 0,
   mHu2IN -> 0,
   Aeij -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   Adij -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   Auij -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mq2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   ml2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   md2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   mu2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   me2Input -> {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
   MassBInput -> 0,
   MassGInput -> 0
};

Options[FSMSSMOpenHandle] = {
    Sequence @@ fsDefaultSettings,
    Sequence @@ fsDefaultSMParameters,
    Sequence @@ fsMSSMDefaultInputParameters
};

FSMSSMOpenHandle[a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSMSSMOpenHandle[a, Sequence @@ s, r];

FSMSSMOpenHandle[OptionsPattern[]] :=
    FSMSSMOpenHandleLib[
        FSMSSMCheckIsNumeric /@ {
            (* spectrum generator settings *)
            OptionValue[precisionGoal],
            OptionValue[maxIterations],
            OptionValue[calculateStandardModelMasses],
            OptionValue[poleMassLoopOrder],
            OptionValue[ewsbLoopOrder],
            OptionValue[betaFunctionLoopOrder],
            OptionValue[thresholdCorrectionsLoopOrder],
            OptionValue[higgs2loopCorrectionAtAs],
            OptionValue[higgs2loopCorrectionAbAs],
            OptionValue[higgs2loopCorrectionAtAt],
            OptionValue[higgs2loopCorrectionAtauAtau],
            OptionValue[forceOutput],
            OptionValue[topPoleQCDCorrections],
            OptionValue[betaZeroThreshold],
            OptionValue[forcePositiveMasses],
            OptionValue[poleMassScale],
            OptionValue[eftPoleMassScale],
            OptionValue[eftMatchingScale],
            OptionValue[eftMatchingLoopOrderUp],
            OptionValue[eftMatchingLoopOrderDown],
            OptionValue[eftHiggsIndex],
            OptionValue[calculateBSMMasses],
            OptionValue[parameterOutputScale],

            (* Standard Model input parameters *)
            OptionValue[alphaEmMZ],
            OptionValue[GF],
            OptionValue[alphaSMZ],
            OptionValue[MZ],
            OptionValue[mbmb],
            OptionValue[Mt],
            OptionValue[Mtau],
            OptionValue[Mv3],
            OptionValue[MW],
            OptionValue[Me],
            OptionValue[Mv1],
            OptionValue[Mm],
            OptionValue[Mv2],
            OptionValue[md2GeV],
            OptionValue[mu2GeV],
            OptionValue[ms2GeV],
            OptionValue[mcmc],
            OptionValue[CKMTheta12],
            OptionValue[CKMTheta13],
            OptionValue[CKMTheta23],
            OptionValue[CKMDelta],
            OptionValue[PMNSTheta12],
            OptionValue[PMNSTheta13],
            OptionValue[PMNSTheta23],
            OptionValue[PMNSDelta],
            OptionValue[PMNSAlpha1],
            OptionValue[PMNSAlpha2],
            OptionValue[alphaEm0],
            OptionValue[Mh]

            (* MSSM input parameters *)
            ,
            OptionValue[TanBeta],
            OptionValue[SignMu],
            OptionValue[MassWBIN],
            OptionValue[Qin],
            OptionValue[QEWSB],
            OptionValue[mHd2IN],
            OptionValue[mHu2IN],
            OptionValue[Aeij][[1,1]],
            OptionValue[Aeij][[1,2]],
            OptionValue[Aeij][[1,3]],
            OptionValue[Aeij][[2,1]],
            OptionValue[Aeij][[2,2]],
            OptionValue[Aeij][[2,3]],
            OptionValue[Aeij][[3,1]],
            OptionValue[Aeij][[3,2]],
            OptionValue[Aeij][[3,3]],
            OptionValue[Adij][[1,1]],
            OptionValue[Adij][[1,2]],
            OptionValue[Adij][[1,3]],
            OptionValue[Adij][[2,1]],
            OptionValue[Adij][[2,2]],
            OptionValue[Adij][[2,3]],
            OptionValue[Adij][[3,1]],
            OptionValue[Adij][[3,2]],
            OptionValue[Adij][[3,3]],
            OptionValue[Auij][[1,1]],
            OptionValue[Auij][[1,2]],
            OptionValue[Auij][[1,3]],
            OptionValue[Auij][[2,1]],
            OptionValue[Auij][[2,2]],
            OptionValue[Auij][[2,3]],
            OptionValue[Auij][[3,1]],
            OptionValue[Auij][[3,2]],
            OptionValue[Auij][[3,3]],
            OptionValue[mq2Input][[1,1]],
            OptionValue[mq2Input][[1,2]],
            OptionValue[mq2Input][[1,3]],
            OptionValue[mq2Input][[2,1]],
            OptionValue[mq2Input][[2,2]],
            OptionValue[mq2Input][[2,3]],
            OptionValue[mq2Input][[3,1]],
            OptionValue[mq2Input][[3,2]],
            OptionValue[mq2Input][[3,3]],
            OptionValue[ml2Input][[1,1]],
            OptionValue[ml2Input][[1,2]],
            OptionValue[ml2Input][[1,3]],
            OptionValue[ml2Input][[2,1]],
            OptionValue[ml2Input][[2,2]],
            OptionValue[ml2Input][[2,3]],
            OptionValue[ml2Input][[3,1]],
            OptionValue[ml2Input][[3,2]],
            OptionValue[ml2Input][[3,3]],
            OptionValue[md2Input][[1,1]],
            OptionValue[md2Input][[1,2]],
            OptionValue[md2Input][[1,3]],
            OptionValue[md2Input][[2,1]],
            OptionValue[md2Input][[2,2]],
            OptionValue[md2Input][[2,3]],
            OptionValue[md2Input][[3,1]],
            OptionValue[md2Input][[3,2]],
            OptionValue[md2Input][[3,3]],
            OptionValue[mu2Input][[1,1]],
            OptionValue[mu2Input][[1,2]],
            OptionValue[mu2Input][[1,3]],
            OptionValue[mu2Input][[2,1]],
            OptionValue[mu2Input][[2,2]],
            OptionValue[mu2Input][[2,3]],
            OptionValue[mu2Input][[3,1]],
            OptionValue[mu2Input][[3,2]],
            OptionValue[mu2Input][[3,3]],
            OptionValue[me2Input][[1,1]],
            OptionValue[me2Input][[1,2]],
            OptionValue[me2Input][[1,3]],
            OptionValue[me2Input][[2,1]],
            OptionValue[me2Input][[2,2]],
            OptionValue[me2Input][[2,3]],
            OptionValue[me2Input][[3,1]],
            OptionValue[me2Input][[3,2]],
            OptionValue[me2Input][[3,3]],
            OptionValue[MassBInput],
            OptionValue[MassGInput]
        }
];

Options[FSMSSMSet] = Options[FSMSSMOpenHandle];

FSMSSMSet[handle_Integer, a___, (fsSettings | fsSMParameters | fsModelParameters) -> s_List, r___] :=
    FSMSSMSet[handle, a, Sequence @@ s, r];

FSMSSMSet[handle_Integer, p:OptionsPattern[]] :=
    FSMSSMSetLib[
        handle,
        ReleaseHold[Hold[{
            (* spectrum generator settings *)
            OptionValue[precisionGoal],
            OptionValue[maxIterations],
            OptionValue[calculateStandardModelMasses],
            OptionValue[poleMassLoopOrder],
            OptionValue[ewsbLoopOrder],
            OptionValue[betaFunctionLoopOrder],
            OptionValue[thresholdCorrectionsLoopOrder],
            OptionValue[higgs2loopCorrectionAtAs],
            OptionValue[higgs2loopCorrectionAbAs],
            OptionValue[higgs2loopCorrectionAtAt],
            OptionValue[higgs2loopCorrectionAtauAtau],
            OptionValue[forceOutput],
            OptionValue[topPoleQCDCorrections],
            OptionValue[betaZeroThreshold],
            OptionValue[forcePositiveMasses],
            OptionValue[poleMassScale],
            OptionValue[eftPoleMassScale],
            OptionValue[eftMatchingScale],
            OptionValue[eftMatchingLoopOrderUp],
            OptionValue[eftMatchingLoopOrderDown],
            OptionValue[eftHiggsIndex],
            OptionValue[calculateBSMMasses],
            OptionValue[parameterOutputScale],

            (* Standard Model input parameters *)
            OptionValue[alphaEmMZ],
            OptionValue[GF],
            OptionValue[alphaSMZ],
            OptionValue[MZ],
            OptionValue[mbmb],
            OptionValue[Mt],
            OptionValue[Mtau],
            OptionValue[Mv3],
            OptionValue[MW],
            OptionValue[Me],
            OptionValue[Mv1],
            OptionValue[Mm],
            OptionValue[Mv2],
            OptionValue[md2GeV],
            OptionValue[mu2GeV],
            OptionValue[ms2GeV],
            OptionValue[mcmc],
            OptionValue[CKMTheta12],
            OptionValue[CKMTheta13],
            OptionValue[CKMTheta23],
            OptionValue[CKMDelta],
            OptionValue[PMNSTheta12],
            OptionValue[PMNSTheta13],
            OptionValue[PMNSTheta23],
            OptionValue[PMNSDelta],
            OptionValue[PMNSAlpha1],
            OptionValue[PMNSAlpha2],
            OptionValue[alphaEm0],
            OptionValue[Mh]

            (* MSSM input parameters *)
            ,
            OptionValue[TanBeta],
            OptionValue[SignMu],
            OptionValue[MassWBIN],
            OptionValue[Qin],
            OptionValue[QEWSB],
            OptionValue[mHd2IN],
            OptionValue[mHu2IN],
            OptionValue[Aeij][[1,1]],
            OptionValue[Aeij][[1,2]],
            OptionValue[Aeij][[1,3]],
            OptionValue[Aeij][[2,1]],
            OptionValue[Aeij][[2,2]],
            OptionValue[Aeij][[2,3]],
            OptionValue[Aeij][[3,1]],
            OptionValue[Aeij][[3,2]],
            OptionValue[Aeij][[3,3]],
            OptionValue[Adij][[1,1]],
            OptionValue[Adij][[1,2]],
            OptionValue[Adij][[1,3]],
            OptionValue[Adij][[2,1]],
            OptionValue[Adij][[2,2]],
            OptionValue[Adij][[2,3]],
            OptionValue[Adij][[3,1]],
            OptionValue[Adij][[3,2]],
            OptionValue[Adij][[3,3]],
            OptionValue[Auij][[1,1]],
            OptionValue[Auij][[1,2]],
            OptionValue[Auij][[1,3]],
            OptionValue[Auij][[2,1]],
            OptionValue[Auij][[2,2]],
            OptionValue[Auij][[2,3]],
            OptionValue[Auij][[3,1]],
            OptionValue[Auij][[3,2]],
            OptionValue[Auij][[3,3]],
            OptionValue[mq2Input][[1,1]],
            OptionValue[mq2Input][[1,2]],
            OptionValue[mq2Input][[1,3]],
            OptionValue[mq2Input][[2,1]],
            OptionValue[mq2Input][[2,2]],
            OptionValue[mq2Input][[2,3]],
            OptionValue[mq2Input][[3,1]],
            OptionValue[mq2Input][[3,2]],
            OptionValue[mq2Input][[3,3]],
            OptionValue[ml2Input][[1,1]],
            OptionValue[ml2Input][[1,2]],
            OptionValue[ml2Input][[1,3]],
            OptionValue[ml2Input][[2,1]],
            OptionValue[ml2Input][[2,2]],
            OptionValue[ml2Input][[2,3]],
            OptionValue[ml2Input][[3,1]],
            OptionValue[ml2Input][[3,2]],
            OptionValue[ml2Input][[3,3]],
            OptionValue[md2Input][[1,1]],
            OptionValue[md2Input][[1,2]],
            OptionValue[md2Input][[1,3]],
            OptionValue[md2Input][[2,1]],
            OptionValue[md2Input][[2,2]],
            OptionValue[md2Input][[2,3]],
            OptionValue[md2Input][[3,1]],
            OptionValue[md2Input][[3,2]],
            OptionValue[md2Input][[3,3]],
            OptionValue[mu2Input][[1,1]],
            OptionValue[mu2Input][[1,2]],
            OptionValue[mu2Input][[1,3]],
            OptionValue[mu2Input][[2,1]],
            OptionValue[mu2Input][[2,2]],
            OptionValue[mu2Input][[2,3]],
            OptionValue[mu2Input][[3,1]],
            OptionValue[mu2Input][[3,2]],
            OptionValue[mu2Input][[3,3]],
            OptionValue[me2Input][[1,1]],
            OptionValue[me2Input][[1,2]],
            OptionValue[me2Input][[1,3]],
            OptionValue[me2Input][[2,1]],
            OptionValue[me2Input][[2,2]],
            OptionValue[me2Input][[2,3]],
            OptionValue[me2Input][[3,1]],
            OptionValue[me2Input][[3,2]],
            OptionValue[me2Input][[3,3]],
            OptionValue[MassBInput],
            OptionValue[MassGInput]
        }] /. HoldPattern[OptionValue[param_]] :> FSMSSMCheckIsNumeric[param] /.
        { p } /.
        FSMSSMGetSettings[handle] /.
        FSMSSMGetSMInputParameters[handle] /.
        FSMSSMGetInputParameters[handle]]];

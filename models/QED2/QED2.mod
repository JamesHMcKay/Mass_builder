(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)
(*                                                                             *)
(*         This file has been automatically generated by FeynRules.            *)
(*                                                                             *)
(* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *)


FR$ModelInformation={
  ModelName->"QED",
  Authors -> {"V. Shtabovenko"},
  Institutions -> {"Technische Universita"t Mu"nchen"},
  Emails -> {"v.shtabovenko@tum.de"},
  Date -> "November 2, 2016"};

FR$ClassesTranslation={};

FR$InteractionOrderPerturbativeExpansion={};

FR$GoldstoneList={};

(*     Declared indices    *)

IndexRange[ Index[Generation] ] = Range[ 3 ]

(*     Declared particles    *)

M$ClassesDescription = {
U[1] == {
    SelfConjugate -> False,
    QuantumNumbers -> {GhostNumber},
    PropagatorLabel -> "uA",
    PropagatorType -> GhostDash,
    PropagatorArrow -> Forward,
    Mass -> 0,
    Indices -> {} },

F[2] == {
    SelfConjugate -> False,
    Indices -> {Index[Generation]},
    QuantumNumbers -> {Q, LeptonNumber},
    PropagatorLabel -> "l",
    PropagatorType -> Straight,
    PropagatorArrow -> Forward*Ghost -> u,
    Mass -> Mlep },

V[1] == {
    SelfConjugate -> True,
    PropagatorLabel -> "\\gamma",
    PropagatorType -> Sine,
    PropagatorArrow -> None,
    Mass -> ma,
    Indices -> {} }
}


(*        Definitions       *)


Mlep[ 1 ] := ME;
Mlep[ 2 ] := MM;
Mlep[ 3 ] := ML;
ma[ ___ ] := ma;


TheLabel[ F[2, {1}] ] := "e";
TheLabel[ F[2, {2}] ] := "mu";
TheLabel[ F[2, {3}] ] := "ta";


(*      Couplings (calculated by FeynRules)      *)

M$CouplingMatrices = {

C[ V[1] , V[1] ] == {{0, (-I)*Z3}, {0, I*Z3}},

C[ -F[2, {e1x2}] , F[2, {e2x2}] , V[1] ] == {{I*EL*IndexDelta[e1x2, e2x2], I*EL*Z1*IndexDelta[e1x2, e2x2]}, {I*EL*IndexDelta[e1x2, e2x2], I*EL*Z1*IndexDelta[e1x2, e2x2]}},

C[ -F[2, {e1x2}] , F[2, {e2x2}] ] == {{0, I*Z2m*IndexDelta[e1x2, e2x2]*Mlep[e1x2]}, {0, I*Z2m*IndexDelta[e1x2, e2x2]*Mlep[e1x2]}, {0, I*Z2z*IndexDelta[e1x2, e2x2]}, {0, I*Z2z*IndexDelta[e1x2, e2x2]}, {0, 0}, {0, 0}}

}


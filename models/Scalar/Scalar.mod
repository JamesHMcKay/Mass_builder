(*
*)
IndexRange[ Index[Generation] ] = Range[1]

M$ClassesDescription = {
  S[1] == {
	SelfConjugate -> True,
	Mass -> ms,
  Indices -> {},
	PropagatorLabel -> "s",
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None }
}
M$CouplingMatrices = {
  (*C[ S[1], S[1]] == -I * {{0,dc}},*)
  C[ S[1], S[1], S[1]] == -I*{{g}},
  C[ S[1], S[1], S[1],S[1]] == -I*{{lambda}}
}
M$LastModelRules = {}

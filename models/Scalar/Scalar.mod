(*
*)
IndexRange[ Index[Generation] ] = Range[1]

M$ClassesDescription = {
  S[1] == {
	SelfConjugate -> True,
	Mass -> Ms,
  Indices -> {},
	PropagatorLabel -> "s",
	PropagatorType -> ScalarDash,
	PropagatorArrow -> None }
}
M$CouplingMatrices = {
  C[ S[1], S[1]] == -I *  { {0, d1Z},{0, d1M} },
  C[ S[1], S[1], S[1]] == -I*{{g,dg}},
  C[ S[1], S[1], S[1],S[1]] == -I*{{lambda,dlambda}}
}
M$LastModelRules = {}

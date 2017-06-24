(* MassBuilder Mathematica package *)
(* Contains functions used during the sorting of amplitdes *)
(* James McKay *)
(* June 2017 *)

BeginPackage["MassBuilder`"]

expandBasisIntegrals::usage =
  "expandBasisIntegrals[Amplitude,Masses] expands the basis integrals into a finite plus a divergent piece.";


makeFiniteAmplitude::usage =
  "makeFiniteAmplitude[amplitude_,order,D_] := extract the terms of specified order in epsilon "

makeFiniteCT::usage =
  "makeFiniteCT[amplitude_,order,D_] := extract the terms of specified order in epsilon after multiplying by 1/epsilon "


MassBuilderEpsilon;
MassBuilderP;

MassBuilderAe[m1_];
MassBuilderBe[m1_,m2_];


Begin["`Private`"]


expandBasisIntegrals[amplitude_, masses_,massesExpand_,A_,B_,J_,K_,T_,V_,F_]:=Module[{amp},amp = amplitude;
amp = expandAIntegrals[amp,masses,massesExpand,A];
amp = expandBIntegrals[amp,masses,B];
amp = expandJKTIntegrals[amp,masses,massesExpand,J,K,T,A];
amp = expandVIntegrals[amp,masses,V,B];
amp];



expandJKTIntegrals[amplitude_, masses_,massesExpand_,J_,K_,T_,A_] := Module[{amp},amp = amplitude;
Do[
   If[
        (massesExpand[[i]]!=0) && (massesExpand[[j]]!=0) && (massesExpand[[k]]!=0),
    amp = amp
    /.
    K[masses[[i]],masses[[j]],masses[[k]],D]  ->
            (
              K[masses[[i]],masses[[j]],masses[[k]],4]
            - (masses[[i]]^2 + masses[[j]]^2 + masses[[k]]^2 )/ (2* MassBuilderEpsilon^2)
            - (masses[[i]]^2 + masses[[j]]^2 + masses[[k]]^2 )/  ( 2*MassBuilderEpsilon)
            + (I*A[masses[[i]],4]+ I*A[masses[[j]],4]+ I*A[masses[[k]],4])/MassBuilderEpsilon
            + MassBuilderAe[masses[[i]]] +  MassBuilderAe[masses[[j]]] + MassBuilderAe[masses[[k]]]
            )];
    amp = amp
    /.
    J[masses[[i]],masses[[j]],masses[[k]],D]  ->
            (
              J[masses[[i]],masses[[j]],masses[[k]],4]
            - (masses[[i]]^2 + masses[[j]]^2 + masses[[k]]^2 )/ (2* MassBuilderEpsilon^2)
            + (-(masses[[i]]^2 +masses[[j]]^2 +masses[[k]]^2)/2 + MassBuilderP^2/4 )/MassBuilderEpsilon
            + ( I*A[masses[[i]],4] + I*A[masses[[j]],4] + I*A[masses[[k]],4] )/MassBuilderEpsilon
            + MassBuilderAe[masses[[i]]] + MassBuilderAe[masses[[j]]] + MassBuilderAe[masses[[k]]]
            );
    If[
        (massesExpand[[i]]!=0),
    amp = amp
    /.
    T[masses[[i]],masses[[j]],masses[[k]],D]  ->
            (
              T[masses[[i]],masses[[j]],masses[[k]],4]
            - 1/(2*MassBuilderEpsilon^2) + 1/(2*MassBuilderEpsilon)
            + ((I*A[masses[[i]],4])/masses[[i]]^2 )/MassBuilderEpsilon
            - ((I*A[masses[[i]],4])/masses[[i]]^2 )
            + MassBuilderAe[masses[[i]]]/ masses[[i]]^2
            )]
            ,
{i, 1,Length[masses]},
{j, 1,Length[masses]},
{k, 1,Length[masses]}];
amp];




expandVIntegrals[amplitude_, masses_,V_,B_] := Module[{amp},amp = amplitude;
Do[
amp = amp
        /.
        V[masses[[i]],masses[[j]],masses[[k]],masses[[l]],D]  ->
        (
          V[masses[[i]],masses[[j]],masses[[k]],masses[[l]],4]
        -1/(2*MassBuilderEpsilon^2) - 1/(2*MassBuilderEpsilon)
        - ( - I * B[masses[[j]] , masses[[l]] ] ) /MassBuilderEpsilon
        - MassBuilderBe[ masses[[j]] , masses[[l]] ]
        ),
{i, 1,Length[masses]},
{j, 1,Length[masses]},
{k, 1,Length[masses]},
{l, 1,Length[masses]}];
amp];


expandBIntegrals[amplitude_, masses_,B_] := Module[{amp},amp = amplitude;
Do[
        amp = amp
        /.
        B[masses[[i]],masses[[j]],D]  ->
        (
          B[masses[[i]],masses[[j]],4]
          + I / MassBuilderEpsilon
          + I * MassBuilderEpsilon * MassBuilderBe[ masses[[i]], masses[[j]] ]
        ),
{i, Length[masses]},
{j, Length[masses]}
];
amp];


expandAIntegrals[amplitude_, masses_,massesExpand_,A_] := Module[{amp},amp = amplitude;
Do[
   If[
        massesExpand[[i]]!=0,
        amp = amp
      /. A[masses[[i]],D]  ->
      (
        A[masses[[i]],4]
      + I * masses[[i]]^2/MassBuilderEpsilon
      - I * MassBuilderEpsilon * MassBuilderAe[masses[[i]]]
      )
      ],
{i, Length[masses]}];
amp];



makeFiniteAmplitude[amplitude_,order_,D_] := Module[{amp,result},amp = amplitude;
                                                         amp = amplitude /. D-> (4 - 2*MassBuilderEpsilon);
                                                         result = Coefficient[amp, MassBuilderEpsilon, order];
                                                         result = result /. MassBuilderEpsilon -> 0;
                                                         result
                                                         ]

makeFiniteCT[amplitude_,order_,D_] := Module[{amp,result},amp = amplitude*(1/MassBuilderEpsilon);
                                                         amp = amp /. D-> (4 - 2*MassBuilderEpsilon);
                                                         result = Coefficient[amp, MassBuilderEpsilon, order];
                                                         result = result /. MassBuilderEpsilon -> 0;
                                                         result
                                                         ]


End[ ]

EndPackage[ ]
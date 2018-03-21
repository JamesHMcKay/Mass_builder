Off[General::spell]

Model`Name = "MDM";
Model`NameLaTechi ="Minimal dark matter";
Model`Authors = "J.McKay";
Model`Date = "2016-08-26";



(*-------------------------------------------*)
(*   Particle Content*)
(*-------------------------------------------*)

(* Gauge Groups *)

Gauge[[1]]={B,   U[1], hypercharge, g1,False};
Gauge[[2]]={WB, SU[2], left,        g2,True};
Gauge[[3]]={G,  SU[3], color,       g3,False};


(* Matter Fields *)

FermionFields[[1]] = {q, 3, {uL, dL},     1/6, 2,  3};  
FermionFields[[2]] = {l, 3, {vL, eL},    -1/2, 2,  1};


FermionFields[[3]] = {d, 3, conj[dR],     1/3, 1, -3};
FermionFields[[4]] = {u, 3, conj[uR],    -2/3, 1, -3};
FermionFields[[5]] = {e, 3, conj[eR],       1, 1,  1};

FermionFields[[6]] = {r, 1,
{{{{ gL, -cL/Sqrt[4]},              {-cL/Sqrt[4], nL/Sqrt[6]}},
{{   -cL/Sqrt[4], nL/Sqrt[6]},      {nL/Sqrt[6],conj[cR]/Sqrt[4]}}},
{{{  -cL/Sqrt[4], nL/Sqrt[6]},      {nL/Sqrt[6], conj[cR]/Sqrt[4]}},
{{   nL/Sqrt[6], conj[cR]/Sqrt[4]}, {conj[cR]/Sqrt[4],conj[gR]}}}}
,0, 5 , 1};


ScalarFields[[1]] =  {H, 1, {Hp, H0},     1/2, 2,  1};



(*----------------------------------------------*)
(*   DEFINITION                                 *)
(*----------------------------------------------*)

NameOfStates={GaugeES, EWSB};

(* ----- Before EWSB ----- *)

DEFINITION[GaugeES][LagrangianInput]= {
	{LagHC, {AddHC->True}},
	{LagNoHC,{AddHC->False}}
};

LagNoHC = mu2 conj[H].H - LamH conj[H].H.conj[H].H;
LagHC =  -(Yd conj[H].d.q + Ye conj[H].e.l + Yu H.u.q)-Yc r.r;
		  		  

(* Gauge Sector *)

DEFINITION[EWSB][GaugeSector] =
{ 
  {{VB,VWB[3]},{VP,VZ},ZZ},
  {{VWB[1],VWB[2]},{VWp,conj[VWp]},ZW}
};     
        
        
          	

(* ----- VEVs ---- *)

DEFINITION[EWSB][VEVs]= 
{    {H0, {v, 1/Sqrt[2]}, {Ah, \[ImaginaryI]/Sqrt[2]},{hh, 1/Sqrt[2]}}     };
 

DEFINITION[EWSB][MatterSector]=   
    {{{{dL}, {conj[dR]}}, {{DL,Vd}, {DR,Ud}}},
     {{{uL}, {conj[uR]}}, {{UL,Vu}, {UR,Uu}}},
     {{{eL}, {conj[eR]}}, {{EL,Ve}, {ER,Ue}}}
     };


(*------------------------------------------------------*)
(* Dirac-Spinors *)
(*------------------------------------------------------*)

DEFINITION[EWSB][Phases] = {
    {Fc, PhaseFc},
    {Fn, PhaseFn},
    {Fg, PhaseFg}
}; 

DEFINITION[EWSB][DiracSpinors]={
 Fd ->{  DL, conj[DR]},
 Fe ->{  EL, conj[ER]},
 Fu ->{  UL, conj[UR]},
 Fv ->{  vL, 0},
 Fc ->{  cL, cR},
 Fg ->{  gL, gR},
 Fn ->{  nL, conj[nL]}};

DEFINITION[EWSB][GaugeES]={
 Fc1 ->{ FcL, 0 },
 Fc2 ->{ 0,  FcR},
 Fd1 ->{  FdL, 0},
 Fd2 ->{  0, FdR},
 Fu1 ->{  Fu1, 0},
 Fu2 ->{  0, Fu2},
 Fe1 ->{  Fe1, 0},
 Fe2 ->{  0, Fe2}};




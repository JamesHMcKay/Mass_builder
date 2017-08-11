(* ::Package:: *)

Quit[]


$LoadTARCER = True;
$LoadFeynArts = True;
<< FeynCalc/FeynCalc.m;
AppendTo[$Path, "/Users/jamesmckay/Documents/Programs/Mass_builder/src/"];
<<MassBuilder.m;


SetOptions[DiracSlash,Dimension->D,FeynCalcInternal->True];SetOptions[DiracTrace,DiracTraceEvaluate->True];null=0;MassBuilderA[mass_,D_]:=TAI[D,0,{{1,mass}}];MassBuilderB[mass1_,mass2_,D_]:=TBI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2}}];MassBuilderJ[mass1_,mass2_,mass3_,D_]:=TJI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2},{1,mass3}}];MassBuilderT[mass1_,mass2_,mass3_,D_]:=TJI[D,Pair[Momentum[p,D],Momentum[p,D]],{{2,mass1},{1,mass2},{1,mass3}}];MassBuilderK[mass1_,mass2_,mass3_,D_]:=TJI[D,0,{{1,mass1},{1,mass2},{1,mass3}}];MassBuilderV[mass1_,mass2_,mass3_,mass4_,D_]:=TVI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2},{1,mass3},{1,mass4}}];MassBuilderF[mass1_,mass2_,mass3_,mass4_,mass5_,D_]:=TFI[D,Pair[Momentum[p,D],Momentum[p,D]],{{1,mass1},{1,mass2},{1,mass3},{1,mass4},{1,mass5}}];FCGV["EL"]=e;FCGV["SW"]=sw;FCGV["CW"]=cw;FCGV["MW"]=mw;FCGV["MZ"]=mz;FCGV["ME"]=me;FCGV["MM"]=mm;FCGV["ML"]=ml;FCGV["MU"]=mu;FCGV["MT"]=mt;FCGV["MD"]=md;FCGV["MS"]=ms;FCGV["MB"]=mb;ZNeu[1,1]=0;ZNeu[1,2]=1;ZNeu[1,3]=0;ZNeu[1,4]=0;ZNeu[2,1]=0;ZNeu[2,2]=0;ZNeu[2,3]=0;ZNeu[2,4]=0;ZNeu[3,1]=0;ZNeu[3,2]=0;ZNeu[3,3]=0;ZNeu[3,4]=0;ZNeu[4,1]=0;ZNeu[4,2]=0;ZNeu[4,3]=0;ZNeu[4,4]=0;UCha[1,1]=1;UCha[1,2]=0;UCha[2,1]=0;UCha[2,2]=0;VCha[1,1]=1;VCha[1,2]=0;VCha[2,1]=0;VCha[2,2]=0;CKM=IndexDelta;MCha[1]=MChi;MCha[2]=MChi;MNeu[1]=MChi;MNeu[2]=MChi;Mh0=mh;SBA=1;


(* ::Subsection:: *)
(* Get the amplitudes that have been computed previously*)
(* Change the path below to the location of your Mass_Builder/models directory *)


path = "/Users/jamesmckay/Documents/Programs/Mass_builder/";
kappa=1/(16\[Pi]^2);
sw=Sin[\[Theta]];
cw=Cos[\[Theta]];

(* Check that counter-term couplings agree with Ibe et al. (2013) *)
g=e/sw;
gammaZ = -e^2(-5/3);
mB0[m1_,m2_]=-I TBI[4, Pair[Momentum[p], Momentum[p]], {{1, m1}, {1, m2}}];
mA0[m_]= -I TAI[4, 0, {{1, m}}];

mB1[m1_,m2_]=-((1/(2 p^2))*(mA0[m2]-mA0[m1]+(p^2+m1^2-m2^2)*mB0[m1,m2]));

mB21[m1_,m2_] = (6*(p^2-m1^2)*mB0[m1,m2]+6*mA0[m1]-6*m1^2+p^2)/(18 p^2);

mB22[m1_,m2_]=-p^2 (mB1[m1,m2]+mB21[m1,m2])-(p^2/4)*mB0[m1,m2]-(1/4)*(m1^2-m2^2)*(mB0[m1,m2]+2*mB1[m1,m2]);

PI[m1_,m2_]=-p^2(mB1[m1,m2]+mB21[m1,m2]);

(* Fermion and lepton sector *) 
V1SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_3_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_4_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_5_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_6_1.mx"}]];
V1SE = V1SE + N SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_7_1.mx"}]];
V1SE = V1SE + N SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_8_1.mx"}]];
V1SE = V1SE + N SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_9_1.mx"}]];
V1SE = V1SE + N SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_10_1.mx"}]];
V1SE = V1SE + N SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_11_1.mx"}]];
V1SE = V1SE + N SelfEnergyFinite*kappa;

V1FermionLepton = makeFiniteAmplitude[V1SE,0,D];
SEIbeFermionLepton =  e^2*(1/(2*Pi^2))*( (-1)^2*(3*PI[mf,mf]) + 3(2/3)^2*(2*PI[mf,mf]+PI[mt,mt]) + 3(1/3)^2*(3*PI[mf,mf]) );
FullSimplify[V1FermionLepton - SEIbeFermionLepton/.Pair[Momentum[p],Momentum[p]]->p^2]



(* Gauge boson sector *) 
V1SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_1_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_2_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_13_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_14_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_15_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_16_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_17_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_18_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

V1Boson = makeFiniteAmplitude[V1SE,0,D];
SEIbeBoson = -  e^2*(3/(4*Pi^2))*( mB22[mw,mw]+ p^2/18  ) -e^2*p^2*(1/(4*Pi^2))*( mB0[mw,mw]);
FullSimplify[V1Boson - SEIbeBoson/.Pair[Momentum[p],Momentum[p]]->p^2]


(* Wino sector *) 
V1SE = 0;
ClearScalarProducts[];

Get[FileNameJoin[{path, "/models/MSSM/output/math_data_V1_12_1.mx"}]];
V1SE = V1SE + SelfEnergyFinite*kappa;

V1Wino = makeFiniteAmplitude[V1SE,0,D];
SEIbeWino = e^2*(1/(2*Pi^2))*( PI[MChi,MChi] );
FullSimplify[V1Wino - SEIbeWino/.Pair[Momentum[p],Momentum[p]]->p^2]


Set @@@ solV1[[1]];


deltaZgammaIbe = - e^2(32/9 *Ng - 5/3)


FullSimplify[dZAA1-deltaZgammaIbe]


SEIbeTotal = FullSimplify[SEIbeWino+SEIbeBoson+SEIbeFermionLepton]


SEIbeDIV = SEIbeTotal/.TBI[4, Pair[Momentum[p], Momentum[p]], {{1, mt}, {1, mt}}]->I/\[Epsilon]/.TBI[4, Pair[Momentum[p], Momentum[p]], {{1, mf}, {1, mf}}]->I/\[Epsilon]\
/.TBI[4, Pair[Momentum[p], Momentum[p]], {{1, MChi}, {1, MChi}}]->I/\[Epsilon]/.TBI[4, Pair[Momentum[p], Momentum[p]], {{1, mw}, {1, mw}}]->I/\[Epsilon]\
/.TAI[4, 0, {{1, mf}}]->I*mf^2/\[Epsilon] /.TAI[4, 0, {{1, MChi}}]->I*MChi^2/\[Epsilon] /.TAI[4, 0, {{1, mw}}]->I*mw^2/\[Epsilon] \
/.TAI[4, 0, {{1, mt}}]->I*mt^2/\[Epsilon];


deltaIbe=FullSimplify[Coefficient[SEIbeDIV,\[Epsilon],-1]]


deltaZgammaIbe = - (e^2/(16*Pi^2))*(32/9 *Ng - 5/3)


FullSimplify[deltaZgammaIbe*kappa-deltaIbe/p^2]

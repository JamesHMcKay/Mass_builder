(* ::Package:: *)

Quit[]


(* :Title: GenerateModelQED													*)

(*
	This software is covered by the GNU General Public License 3.
	Copyright (C) 1990-2016 Rolf Mertig
	Copyright (C) 1997-2016 Frederik Orellana
	Copyright (C) 2014-2016 Vladyslav Shtabovenko
*)

(* :Summary:  Generates FeynArts model for QED 								*)

(* ------------------------------------------------------------------------ *)



(* ::Section:: *)
(*QED model for FeynArts*)


(* ::Subsection:: *)
(*Load FeynRules*)


FR$Parallel=False;
$FeynRulesPath=SetDirectory["/Users/jamesmckay/Documents/Programs/FeynRules"];
<<FeynRules`;


(* ::Subsection:: *)
(*Load FeynRules model*)


frModelPath=FileNameJoin[{"/Users/jamesmckay/Documents/Programs/Mass_builder/FeynRules/","QED","QED2.fr"}];
LoadModel[frModelPath];




(* ::Subsection:: *)
(*Create FeynArts model*)


FR$Loop=True;
SetDirectory[FileNameJoin[{"/Users/jamesmckay/Documents","Programs","Mass_builder","models","QED2"}]];
WriteFeynArtsOutput[LQED,Output->"QED",CouplingRename->False];










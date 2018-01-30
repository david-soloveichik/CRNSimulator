(* ::Package:: *)

(* ::Text:: *)
(*Chemical Reaction Network (CRN) Simulator package is developed by David Soloveichik. Copyright 2009-2013. *)
(*http://www.dna.caltech.edu/~davids/*)


Needs["CRNSimulator`"]
BeginPackage["CRNSimulatorExtensions`", {"CRNSimulator`"}];


constflux::usage="Represents a constant flux reaction approximation.";
SimulateRxnsysWithSchedule::usage="";
SteppizeSchedule::usage="";


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


(* ::Subsection:: *)
(*Constant Rate Reactions*)


constflux[r_,p_,k_]:=Module[{f},Seq[rxn[Evaluate[r+f],p,10^($MachinePrecision)],rxn[1,f,k]]]


(* ::Subsection:: *)
(*Controlling Simulations*)


(*Given a solution, generates a new reaction system with the initial concentrations given by the final concentrations of the solution.*)
rsysForRestart[rsys_,sol_,tmax_]:=
Join[
rsys/.conc[___]->Seq[],
conc[#,(#[tmax]/.sol)]&/@SpeciesInRxnsys[rsys]]
(*For formal functions of one argument, whose definitions are valid starting at 0 ending at end1 and end2.  Make new function which uses the second one when first one ends:*)
joinFunctions[f1_,f2_,end1_,end2_]:=
Function[Evaluate[Which@@{0<=#<end1, f1[#],end1<=#<end1+end2, f2[#-end1]}]]
(*Main function for simulating a reaction system using a schedule of addition/removal of species, rather than just starting in one initial configuration.  (In rxnsys, x1 must not have initial conc):*)
SimulateRxnsysWithSchedule[rxnsys_,schedule_,x1_,opts:OptionsPattern[NDSolve]]:=
SimulateRxnsysWithSchedule[rxnsys,schedule,{x1},opts]
SimulateRxnsysWithSchedule::undershoot="Scheduled decrease (schedule step `1`: `2`) resulted in negative concentration. Going to 0.";
SimulateRxnsysWithSchedule[rxnsys_,schedule_,xs_List,opts:OptionsPattern[NDSolve]]:=
Module[
{n=Length[schedule],rsys,lastsol,sol},

rsys=Join[rxnsys, MapIndexed[conc[#1,schedule[[1,1+#2[[1]]]]]&,xs]];
lastsol=SimulateRxnsys[rsys,schedule[[1,1]],opts];
sol=lastsol;

Do[
rsys=
rsysForRestart[rsys,lastsol,schedule[[i-1,1]]]/.
MapIndexed[
conc[#1,s_]:>
(If[s+schedule[[i,1+#2[[1]]]]<0,
Message[SimulateRxnsysWithSchedule::undershoot,i,schedule[[i,1+#2[[1]]]]]];
conc[#1,Max[0,s+schedule[[i,1+#2[[1]]]]]])&,
xs];
lastsol=SimulateRxnsys[rsys,schedule[[i,1]],opts];
sol=MapThread[#1/.(s_->f1_):>(s->joinFunctions[f1,#2/.(_->f2_):>f2,Total[schedule[[;;i-1,1]]],Total[schedule[[;;i,1]]]])&,
{sol,lastsol}],

{i,2,n}];

sol
]

(*Since some species are buffered, we may not be able to instantaneously decrease concentration because this may result in going negative.  We do 1/3 of the distance, then 1/3 of the remaining distance, etc (here 3 is the stepfactor).  Rescaled so that after all iterations we cover the whole distance.*)
SteppizeSchedule[schedule_,iterations_,stepfactor_:3,stepduration_:0.01]:=
Module[{sum=1-((stepfactor-1)/(stepfactor))^(iterations)},
Map[
Seq@@Table[Prepend[N[Rest[#]]* (stepfactor-1)^(i-1)(stepfactor)^(-i) / sum,If[i==iterations,First[#],stepduration]],{i,1,iterations}]&,
schedule]]


(* ::Subsection:: *)
(*Epilog*)


End[];
EndPackage[];

(* ::Package:: *)

(* ::Text:: *)
(*Chemical Reaction Network (CRN) Simulator package is developed by David Soloveichik. Copyright 2009-2015. *)
(*http://users.ece.utexas.edu/~soloveichik/crnsimulator.html*)


Needs["CRNSimulator`"];
BeginPackage["CRNSimulatorSSA`", {"CRNSimulator`"}];


(* ::Section:: *)
(*Public interface specification*)


SimulateRxnsysSSA::usage="SimulateRxnsysSSA[rxnsys,endtime] runs Gillespie's SSA simulation of the reaction \
system rxnsys for time 0 to endtime or until no further reactions are possible. Set endtime to \[Infinity] to continue \
until no further reactions are possible. In rxnsys, reactions are specified by rxn statements (arbitrary \
order reactions are allowed: eg rxn[a+a+b,c+c,1], including zero-order reactions), and initial molecular counts \
are specified by conc statements. \
If no initial count is set for a species, its initial count is set to 0. SimulateRxnsysSSA returns a \
TimeSeries object of reaction firing times and the corresponding molecular counts. The \
propensity of reaction a\[Rule]\[Ellipsis] is k*a, the propensity of reaction a+b\[Rule]\[Ellipsis] is k*a*b, and the propensity of \
reaction a+a\[Rule]\[Ellipsis] is k*a*(a-1), etc.";


(* ::Section:: *)
(*Private*)


Begin["`Private`"];


(* Set Compiler options *)
On["CompilerWarnings"];
System`SetSystemOptions["CompileOptions"->"CompileReportExternal"->True];
On[Compile::noinfo];
SetOptions[Compile,CompilationOptions->{"InlineExternalDefinitions" -> True,"InlineCompiledFunctions"->False},
CompilationTarget->"C",Parallelization->True,RuntimeOptions->"Quality"];


(* Produces reactions in srxn[[reaction index],{x1,x2},{x3},k] format ("structured reaction system"). Removes any other statements. *)
rxnsysToSrxnsys[rsys_] := 
	Module[{i = 1}, 
		Cases[
			rsys, 
			rxn[rs_, ps_, k_] :> 
				srxn[i++, 
					{Unevaluated[rs] /. {Plus -> Seq, c_Integer*s_ :> Seq @@ Table[s, {c}]}}, {Unevaluated[ps] /. {Plus -> Seq, c_Integer*s_ :> Seq @@ Table[s, {c}]}}, 
					k]]]


(* define Gillespie kinetics propensities *)

(* Logically this is: #X!/(#X-mult)! where mult is the multiplicity. But we optimize numerics. 
Note: "if statement" is not logically necessary, but avoids unnecessary multiplication *)
reactantPropensityTerm[s_,mult_]:=If[s>=mult,Times@@(s-Range[0,mult-1]),0]

(* General case definition: Can handle arbitrary arity and multiplicity *)
reactantListToPropensity[rl_]:=Times@@Cases[Tally[rl],{s_,mult_}:>reactantPropensityTerm[s,mult]]


rxnsysToPropensities[rsys_]:=Cases[rxnsysToSrxnsys[rsys],srxn[_,rl_,_,k_]:>k*reactantListToPropensity[rl]]


rxnsysToStateChangeVectors[rsys_] := 
	Module[
		{spcs = Sort[SpeciesInRxnsys[rsys]], 
			srxnsys = rxnsysToSrxnsys[rsys], 
			symbolicNetChanges}, 
		symbolicNetChanges = Cases[srxnsys, srxn[_, r_, p_, _] :> Plus @@ p - Plus @@ r];
		Outer[Coefficient[#1, #2] &, symbolicNetChanges, spcs]
	]


SimulateRxnsysSSA[rsys_, tlimitarg_, maxNumStepsarg_: 10^6] :=
 
 Module[{concs, spcs, x0arg, stateChangeVectors, SSA, w, i, expRvPool,
    unifRvPool, trj, times, data, timesrun, datarun},
  
  spcs = SpeciesInRxnsys[rsys];
  
  (* Vector of initial counts from parsing conc statements. 
  Multiple conc for same species are summed. *)
  
  x0arg = Plus @@ Cases[rsys, conc[#, c_] :> c] & /@ spcs;
  
  stateChangeVectors = rxnsysToStateChangeVectors[rsys];
  
  SSA = With[{computePropensities = 
      Hold[Compile][{{w, _Integer, 1}}, 
        rxnsysToPropensities[rsys] /. 
         MapIndexed[#1 -> (Hold[w[[i]]] /. (i -> (First[#2]))) &, 
          spcs]] // ReleaseHold},
    
    Compile[{{stateChangeVectors, _Integer, 2}, {x0, _Integer, 
       1}, {t0, _Real}, {tlimit, _Real}, {expRvPool, _Real, 
       1}, {unifRvPool, _Real, 1}},
     
     Module[{x = x0, propensities, sumpropensities, t = t0, 
       rvPointer = 1, rxnNum = Length[stateChangeVectors], 
       spcsNum = Length[First[stateChangeVectors]], nextRxn, 
       maxNumSteps = Length[expRvPool], deltat, output, step = 1, s, 
       p = 0.0, mark},
      
      output = Table[0.0, {maxNumSteps}, {spcsNum + 1}];  (* 
      the first dimension is time *)
      
      
      While[
       If[tlimit == 0, True, t <= tlimit] && 
        If[maxNumSteps == 0, True, step <= maxNumSteps], (* 
       use tlimit or maxNumSteps only if nonzero *)
       
       propensities = computePropensities[x];
       sumpropensities = Total[propensities];
       If[sumpropensities == 0.0, Break[]];  (* 
       no reaction is possible *)
       
       (* compute deltat and nextRxn*)
       
       deltat = (1/sumpropensities)*expRvPool[[rvPointer]];
       mark = sumpropensities*unifRvPool[[rvPointer]];
       rvPointer++;
       For[p = 0.0; nextRxn = 1, nextRxn <= rxnNum, nextRxn++,
        p += propensities[[nextRxn]];
        If[p >= mark, Break[]]
        ];
       
       (* implement reaction *)
       
       x += stateChangeVectors[[nextRxn]];
       t += deltat;
       
       output[[step, 1]] = t;
       For[s = 1, s <= spcsNum, s++,
        output[[step, 1 + s]] = x[[s]];
        ];
       
       step++;
       ];
      output[[;; step - 1]]
      ]
     ]
    ];
  
  
  (* initial state as first point *)
  times = {0.0};
  data = {x0arg};
  
  While[True,
   expRvPool = 
    RandomVariate[ExponentialDistribution[1], maxNumStepsarg];
   unifRvPool = RandomVariate[UniformDistribution[], maxNumStepsarg];
   trj = SSA[stateChangeVectors, Last[data], Last[times], 
     tlimitarg /. \[Infinity] -> 0, expRvPool, unifRvPool];
   timesrun = trj[[All, 1]];
   datarun = Round[trj[[All, 2 ;;]]];
   times = Join[times, timesrun];
   data = Join[data, datarun];
   If[Length[timesrun] < maxNumStepsarg, Break[]];
   ];
  
  TimeSeries[data,{times},ResamplingMethod->{"Interpolation",InterpolationOrder->0}]
  
]



(* Reset Compile options so Plot, etc works ok *)
SetOptions[Compile,CompilationOptions->All];


End[];
EndPackage[];

(* ::Package:: *)

(* ::Package:: *)
(* :Title: FeynFlow *)
(* :Context: FeynFlow` *)
(* :Author: Reefat, A. Aleksejevs, and S. Barkanova *)


(* :Summary: FeynFlow: Automated Feynman Rule Derivation and Model File Generation. 
             If used, please cite: Reefat, A. Aleksejevs, and S. Barkanova, 
             FeynFlow: Automated Feynman Rule Derivation and Model File Generation, 
             Zenodo DOI: https://doi.org/10.5281/zenodo.22285890 (2026). *)
(* :Package Version: 1.0 *)
(* :Mathematica Version: 12.2 *)

BeginPackage["FeynFlow`", {"FeynCalc`"}]

GetCounterTermLagrangian::usage = 
  "GetCounterTermLagrangian[inputLagrangian, scalarList] generates the counter term Lagrangian using multiplicative renormalization.";

FAHReplace::usage = 
  "FAHReplace[exp, changes] replaces delta Z symbols with output-friendly format.";

FAHFeynmanRules::usage = 
  "FAHFeynmanRules[Lagrangian, fieldCombination] generates Feynman rules for a given Lagrangian.";

FAHCTFeynmanRules::usage = 
  "FAHCTFeynmanRules[Lagrangian, scalars, fieldCombination] generates Feynman rules for counter term Lagrangian.";

GetCouplingVector::usage = 
  "GetCouplingVector[treeLevelexpr, loopLevelexpr, lorentzStructure] extracts coupling vectors from tree and loop level expressions.";

WriteCV::usage = 
  "WriteCV[fieldCombo, cv, fileName] writes coupling vectors to a FeynArts compatible .mod file.";

GetPropagator::usage = 
  "GetPropagator[exp, lorStructure, repRules] calculates the propagator from given expression.";

GenPropInternal::usage = 
  "GenPropInternal[texp, field, gaugeSym, lorRep] generates internal propagator in FeynArts format.";

WriteGenProp::usage = 
  "WriteGenProp[expr, fileName] writes generic propagators to file.";

WriteGENOutput::usage = 
  "WriteGENOutput[texp, lorRep, fileName] writes Lorentz structures in FeynArts output format.";

Begin["`Private`"]

(* Ensure FeynCalc context is accessible *)
If[Not[MemberQ[$ContextPath, "FeynCalc`"]],
  PrependTo[$ContextPath, "FeynCalc`"]
];

If[Not[MemberQ[$Packages, "FeynCalc`"]],
  Quiet[Needs["FeynCalc`"]]
];

(* Helper function for complete simplification *)
FullFeynCalcSimplify[expr_] := expr // 
  DiracSimplify // 
  Contract // 
  DiracSimplify // 
  SMPToSymbol//
  FCE // 
  FCI//
  DiracGamma//
  Momentum//
  DiracSimplify // 
  FCReplaceAll//
  Simplify // 
  Collect2[#, {DiracGamma, Pair}] &;

(* ========================================================== *)
(*Operator to get the counter term Lagrangian part*)
GetCounterTermLagrangian[inputLagrangian_,scalarList_List]:=Module[
(***
  arguments are:
	-inputLargrangian: using FeynCalc Opertors 
   -scalarList: expects a List with symbols
return: an expression for CT lagrangian in FeynCalc language. 	
***)
{replacementRules, scalarReplacements},

(*Adding all the replacement rules for the fields*)
replacementRules={
QuantumField[field1_]:> QuantumField[field1]+(1/2)*\[Delta]*SMPToSymbol[Superscript[Subscript[\[Delta]z,2],field1]]*QuantumField[field1],QuantumField[field1_,{lor1_}]:> QuantumField[field1,{lor1}]+(1/2)*\[Delta]*SMPToSymbol[Superscript[Subscript[\[Delta]z,2],field1]]*QuantumField[field1,{lor1}],QuantumField[FCPartialD[lor1_],field1_,{lor2_}]:> QuantumField[FCPartialD[lor1],field1,{lor2}]+(1/2)*\[Delta]*SMPToSymbol[Superscript[Subscript[\[Delta]z,2],field1]]*QuantumField[FCPartialD[lor1],field1,{lor2}],QuantumField[FCPartialD[lor1_],field_]:> QuantumField[FCPartialD[lor1],field]+(1/2)*\[Delta]*SMPToSymbol[Superscript[Subscript[\[Delta]z,2],field]]*QuantumField[FCPartialD[lor1],field]
};

(*Generating all the replacement rules for the scalars given by user
and adding them to the replacementRules array.
*)
scalarReplacements={#:> #+#*\[Delta]*SMPToSymbol[Superscript[Subscript[\[Delta]z,1],#]]&/@scalarList}//Flatten;
replacementRules=Join[replacementRules,scalarReplacements]//Flatten;
(*Print[replacementRules];*)
(*Exapnding the given Lagrangian and applying all the replacement rules
in the Lagrangian to get the Counter term Lagrangian using the 
multiplicative renormalization.
*)

Return[
((Expand[inputLagrangian/.replacementRules]//.\[Delta]^2->0//.\[Delta]^3->0//.\[Delta]^4->0)-Expand[inputLagrangian])//.\[Delta]->1
]

]
(*Operator to change the look of the output that will be easier to process into text files*)
FAHReplace[exp_,changes_Association]:=Module[
{res=exp},
(*SetAttributes[symbolNameJoin,HoldAll];
symbolNameJoin[symbols__Symbol]:=Symbol@Apply[StringJoin,Map[Function[s,ToString[Unevaluated[s]],HoldFirst],Hold[symbols]]];*)

res=res//.SMPToSymbol[Superscript[Subscript[\[Delta]z,1],a_]]:> OutputForm[StringJoin["dZ1",ToString[a]]];

For[i=1,i<=Length[changes],i++,
(*Print[Keys[changes][[i]],"\n",Values[changes][[i]]];*)
res = res//.SMPToSymbol[Superscript[Subscript[\[Delta]z,2],Keys[changes][[i]]]]:> Symbol["dZ2"<>ToString[Values[changes][[i]]]]
];
Return[res]
]
(*Operator to get the Feynman Rules for a lagrangian*)
FAHFeynmanRules[Lagrangian_,fieldCombination_]:=Module[{lag=Lagrangian,possibleCombo=fieldCombination,temp={}},
(***
	Generate the Feynman rules for a given Lagrangian and all possible field combinations.
***)
For[i=1,i<=Length[possibleCombo],i++,					AppendTo[temp,Simplify[Contract[DiracSimplify[Expand[FunctionalD[lag,possibleCombo[[i]]]]
							//.QuantumField[a_,{b_}]:>0 //.QuantumField[FCPartialD[a_],b_,{c_}]:>0 ]]]
				];
		];
	Return[temp]
		];
(*Operator to get the Feynman Rules for the Counter term lagrangian*)
FAHCTFeynmanRules[Lagrangian_,scalars_,fieldCombination_]:=Module[{lag,possibleCombo=fieldCombination,tScalar=scalars,temp={}},
(***
	Generate the Feynman rules for a given Lagrangian and all possible field combinations.
***)
lag=GetCounterTermLagrangian[Lagrangian,tScalar]//Simplify;
	For[i=1,i<=Length[possibleCombo],i++,			AppendTo[temp,Simplify[Contract[DiracSimplify[Expand[FunctionalD[lag,possibleCombo[[i]]]]
							//.QuantumField[a_,{b_}]:>0 //.QuantumField[FCPartialD[a_],b_,{c_}]:>0 ]]]
				];
		];
	Return[temp]
		];
GetCouplingVector[treeLevelexpr_,loopLevelexpr_,lorentzStructure_]:=Module[
{convertedLorStr={},
revertedLorStr={},
wexprTree,wexprLoop,res,f},
(*
Changing the user given Lorentz structure into coefficients
*)
For[i=1,i<=Length[lorentzStructure],i++,
AppendTo[convertedLorStr,FCI[lorentzStructure[[i]]]-> "coef"<>ToString[i]];
];
(*Print[convertedLorStr];*)
(*
Creating another array to get the inversion for final replacements.
*)
For[i=1,i<=Length[convertedLorStr],i++,
AppendTo[revertedLorStr,convertedLorStr[[i]][[2]]->convertedLorStr[[i]][[1]]];
];
(*Print[revertedLorStr];*)
(*
Replacing all the 'coef'-s into the tree level expression and
loop level expression.
*)
wexprTree= ReplaceAll[Expand[treeLevelexpr],convertedLorStr];
wexprLoop= FCReplaceAll[Expand[Calc[loopLevelexpr]*(ChiralityProjector[+1]+ChiralityProjector[-1])],convertedLorStr];
(*Print[wexprTree];
Print[wexprLoop];*)

(*
Creating an array to hold the results. NOTE: 2 is because we have
a tree level and 1-loop in the coupling vector.
*)
res=Array[f,{Length[lorentzStructure],2}];

For[i=1,i<=Length[lorentzStructure],i++,
res[[i]][[1]]=Coefficient[wexprTree,"coef"<>ToString[i]]//.1->0//.-1->0//DiracSimplify//Simplify;
res[[i]][[2]]=Coefficient[wexprLoop,"coef"<>ToString[i]]//.DiracGamma[6]-> (1-DiracGamma[7])//DiracSimplify//Simplify;
];

Return[ReplaceAll[res,revertedLorStr]]
]
(*
WriteCV is an operator that will take in field combination list, coupling vector array and a file name. This will return a *.mod file
where the coupling vectors will be written in FeynArts compatible format.
*)
WriteCV[fieldCombo_List,cv_List,fileName_String]:=Module[
{t1={},t2={},str},
(*
Seperate all possible field combinations and store them in t1.
*)
For[i=1,i<=Length[fieldCombo],i++,
AppendTo[t1,fieldCombo[[i]]]
];

(*
Create the lhs of the Coupling vectors and add them to t2.
*)
For[i=1,i<=Length[fieldCombo],i++,
AppendTo[t2,"C["<>StringJoin[ToString/@Riffle[t1[[i]],","]]<>"]=="]
];

(*
Writing the output in a *.mod file.
*)
str=OpenWrite[FileNameJoin[{Directory[],fileName}]];

Write[str,OutputForm["M$CouplingMatrices = { \n"]];
For[i=1,i<=Length[t1],i++,
If[i<Length[t1],Write[str,OutputForm["C["<>StringJoin[ToString/@Riffle[t1[[i]],","]]<>"]=="],cv[[i]],OutputForm[", \n \n"]],
Write[str,OutputForm["C["<>StringJoin[ToString/@Riffle[t1[[i]],","]]<>"]=="],cv[[i]],OutputForm["\n \n"]]
];
];
Write[str,OutputForm["}"]];
Close[str]
]
(*Propagator Part*)
GetPropagator[exp_,lorStructure_List,repRules_List]:=Module[
{
coefList={},expWithRep=exp/.repRules,listOfEqn={},
quadFormRep,dabPbd,newLorStructure,temResult={},solution,
rrDen
},

(*Create a list of coefficients a1,a2,...,an*)
For[i=1,i<=Length[lorStructure],i++,
AppendTo[coefList,"a"<>ToString[i]];
];

(*Creating a quadratic form with the replacement rules*)
quadFormRep=Total[coefList*lorStructure];

(*Dot Product between given expression and the quadFormRepl*)
dabPbd=((quadFormRep/.repRules)*exp//DiracSimplify//Contract)/.repRules//Simplify;

(*Change of the given Lorentz structure to match the lorentz indices*)
newLorStructure=(lorStructure/.repRules);

(*Creating a List of Equations to be solved later*)
For[i=1,i<=Length[lorStructure],i++,
If[i==1,AppendTo[listOfEqn, Coefficient[FCI[dabPbd],FCI[newLorStructure[[i]]]]==1],
AppendTo[listOfEqn, Coefficient[FCI[dabPbd],FCI[newLorStructure[[i]]]]==0]]
];

(*Solving for the coefficients a1, a2, ... an*)
solution=Solve2[listOfEqn,coefList];

(*Replacing the solutions to get the inverse of
the quadrartic form.
*)
For[i=1,i<=Length[solution],i++,
AppendTo[
temResult,solution[[i]][[2]]*newLorStructure[[i]]
]
];
(*We do not want any Gamma Matrices/ Gamma Slash in the denominator
so we have to introduce the following Denominator Replacement*)
rrDen={(1/(m_ - DiracGamma[Momentum[k_]])):>(m + DiracGamma[Momentum[k]])/(Superscript[m,2] - Pair[Momentum[k], Momentum[k]])};
(*Summing up the list to get the propagator*)
Return[Total[temResult]//.rrDen]
]
GenPropInternal[texp_,field_,gaugeSym_:Null,lorRep_List]:=Module[{exp=(texp/.lorRep),replacementRules,denominatorPart,numeratorPart,m},

replacementRulesAll={gaugeSym:>OutputForm["Sqrt[GaugeXi["<>ToString[field]<>"]]"],((num__)/(m_^2-Pair[Momentum[k_],Momentum[k_]])):>num OutputForm["PropagatorDenomiator["<>ToString[k]<>",Mass["<>ToString[field]<>"]]"],
((num__)/(Superscript[m_,2]-Pair[Momentum[k_],Momentum[k_]])):>num OutputForm["PropagatorDenomiator["<>ToString[k]<>",Mass["<>ToString[field]<>"]]"],
((num__)/(lam_*m_^2-Pair[Momentum[k_],Momentum[k_]])):>num OutputForm["PropagatorDenomiator["<>ToString[k]<>","<>ToString[Sqrt[lam]]<>"Mass["<>ToString[field]<>"]"],(GA[a_]GA[b_]):>NonCommutative[GA[a],GA[b]],(GA[a_]FV[b__]):>NonCommutative[GA[a],FV[b]],(*(FV[a__]FV[b__]):>NonCommutative[FV[a],FV[b]],*)
(DiracGamma[Momentum[a_]]):>OutputForm["DiracSlash["<>ToString[a]<>"]"],
(GS[a_]):>OutputForm["DiracSlash["<>ToString[a]<>"]"],
MT[a_,b_]:>OutputForm["MetricTensor["<>ToString[a]<>","<>ToString[b]<>"]"],FV[a_,b_]:>OutputForm["FourVector["<>ToString[a]<>","<>ToString[b]<>"]"],GA[a_]:>OutputForm["DiracMatrix["<>ToString[a]<>"]"],
(num__/Superscript[m_,2]):>num*OutputForm["Mass["<>ToString[field]<>"]^2"],
(num__/m_^2):>num*OutputForm["Mass["<>ToString[field]<>"]^2"],(num__/m_):>num*OutputForm["Mass["<>ToString[field]<>"]"]};

replacementRulesOnce={(DiracGamma[Momentum[k__]]+c__):>NonCommutative[DiracGamma[Momentum[k]]+c]};


Return[exp/.replacementRulesOnce//.replacementRulesAll]

]
WriteGenProp[expr_List,fileName_]:=Module[
{str},
str=OpenWrite[FileNameJoin[{Directory[],fileName}]];
Write[str,OutputForm["M$GenericPropagators = {"]];
For[i=1,i<=Length[expr],i++,
Write[str, OutputForm["AnalyticalPropagator[Internal][...]=="] ];
If[
i==Length[expr],Write[str,expr[[i]]],
Write[str,expr[[i]],OutputForm[", \n "]]
]
];
Write[str,OutputForm["}"]];
Close[str]
]
(*Lorentz Structure part of the couplings into Output*)
WriteGENOutput[texp_,lorRep_List,fileName_String:Null]:=Module[{exp=(texp/.lorRep),replacementRulesAll,replacementRules,denominatorPart,numeratorPart,m,str,expr},

replacementRulesAll={
(GA[a_]*GA[b_]):>OutputForm["NonCommutative[DiracGamma["<>ToString[a]<>"],DiracGamma["<>ToString[b]<>"]"],
(GA[a_] . GS[p_]*GA[5]):>OutputForm["NonCommutative[DiracMatrix["<>ToString[a]<>"],DiracSlash["<>ToString[p]<>"],(ChiralityProjector[+1]-ChiralityProjector[-1])]"],
(GS[p_] . GA[a_]*GA[5]):>OutputForm["NonCommutative[DiracSlash["<>ToString[p]<>"],(ChiralityProjector[+1]-ChiralityProjector[-1]),DiracMatrix["<>ToString[a]<>"]]"],
GS[a_]:>OutputForm["NonCommutative[DiracSlash["<>ToString[a]<>"]]"],
MT[a_,b_]:>OutputForm["MetricTensor["<>ToString[a]<>","<>ToString[b]<>"]"],FV[a_,b_]:>OutputForm["FourVector["<>ToString[a]<>","<>ToString[b]<>"]"],GA[a_]:>OutputForm["NonCommutative[DiracMatrix["<>ToString[a]<>"]]"],
SP[a_,b_]:>OutputForm["ScalarProduct["<>ToString[a]<>","<>ToString[b]<>"]"],
(DiracGamma[6]*GS[p_]):>OutputForm["NonCommutative[ChiralityProjector[+1],DiracSlash["<>ToString[p]<>"]]"],
(DiracGamma[7]*GS[p_]):>OutputForm["NonCommutative[ChiralityProjector[-1],DiracSlash["<>ToString[p]<>"]]"],
DiracGamma[6]-> OutputForm["NonCommutative[ChiralityProjector[+1]]"],
DiracGamma[7]-> OutputForm["NonCommutative[ChiralityProjector[-1]]"]};

replacementRulesOnce={(DiracGamma[Momentum[k__]]+c__):>NonCommutative[DiracGamma[Momentum[k]]+c]};


expr =exp/.replacementRulesOnce//.replacementRulesAll;
(*Print[expr]*)
(*The following is used to generate the text file*)
If[fileName===Null,
Return[expr],
Print["Writing "<>ToString[fileName]<>"file..."];
str=OpenAppend[FileNameJoin[{Directory[],fileName}],PageWidth->Infinity];
Write[str,OutputForm["AnalyticalCoupling[...]==G[...][...].\n"],expr];
Write[str,OutputForm["\n "]];
Close[str]
];
]
(* ========================================================== *)
End[]

EndPackage[]

Print["FeynFlow loaded successfully."];
Print["If you use this package, please cite:"];
Print["Reefat, A. Aleksejevs, and S. Barkanova,"];
Print["FeynFlow: Automated Feynman Rule Derivation and Model File Generation,"];
Print["Zenodo DOI: https://doi.org/10.5281/zenodo.22285890 (2026)."];

#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Hofstadter Model*)


(* ::Text:: *)
(*Working in Landau gauge that Subscript[A, 1]=-B Subscript[x, 2] and Subscript[A, 2]=0, and assuming that the flux piercing through each plaquette is a fraction of the quantum flux: B\[CenterDot]a^2\[Congruent](p/q)Subscript[\[Phi], 0], then the non-interacting TB Hamiltonian *)


(* ::DisplayFormula:: *)
(*H=-t\!\( *)
(*\*SubscriptBox[\(\[Sum]\), \(r; i = 1, 2\)]\( *)
(*\*SuperscriptBox[\(E\), \(I\ *)
(*\*SubscriptBox[\(A\), \(i\)] \((x)\)\)] \**)
(*TemplateBox[{"x"},*)
(*"Ket"] \**)
(*TemplateBox[{RowBox[{"x", "+", SubscriptBox["x", "i"]}]},*)
(*"Bra"]\)\)+h . c .*)


(* ::Text:: *)
(*gives an iterative Schrodinger E.q. for the on-site Wannier functions:*)


(* ::DisplayFormulaNumbered:: *)
(*-t(E^(-I 2\[Pi](p/q)Subscript[x, 2]) \[Psi](Subscript[x, 1]+1,Subscript[x, 2])+E^(I 2\[Pi](p/q)Subscript[x, 2]) \[Psi](Subscript[x, 2]-1,Subscript[x, 1])+\[Psi](Subscript[x, 1],Subscript[x, 2]+1)+\[Psi](Subscript[x, 1],Subscript[x, 2]-1))=E \[Psi] (Subscript[x, 1],Subscript[x, 2]) .*)


(* ::Text:: *)
(*Clearly we have a enlarged magnetic translation symmetry: (Subscript[x, 1],Subscript[x, 2])->(Subscript[x, 1]+q,Subscript[x, 2]), and (Subscript[x, 1],Subscript[x, 2])->(Subscript[x, 1],Subscript[x, 2]+1), from which we can define the magnetic translation operator Subscript[T, 1] and Subscript[T, 2], such that*)


(* ::DisplayFormulaNumbered:: *)
(*[Subscript[T, 1]^q,Subscript[T, 2]]=[Subscript[T, 1]^q,H]=[Subscript[T, 2],H]=0.*)


(* ::Text:: *)
(*Now we can perform the Fourier transformation to the magnetic Brillouin zone (since the real space translation symmetry is enlarged along Subscript[x, 1]-direction, the corresponding BZ is then shrink to 2\[Pi]/q along the Subscript[k, 1]-direction. And since momentum eigenstates Ket[Subscript[k, 1],Subscript[k, 2]] are eigenstates of \!\(\*SubsuperscriptBox[\(T\), \(1\), \(q\)]\) and Subscript[T, 2], we must have q Landau bands for each Wannier state:*)


(* ::DisplayFormula:: *)
(*\[Psi](Subscript[x, 1],Subscript[x, 2])=1/q \!\( *)
(*\*SubsuperscriptBox[\(\[Sum]\), \(r = 1\), \(q\)]\( *)
(*\*SubsuperscriptBox[\(\[Integral]\), \(-\[Pi]\), \(\[Pi]\)]*)
(*\*FractionBox[\(\[DifferentialD]*)
(*\*SubscriptBox[\(k\), \(1\)]\), \(2  \[Pi]\)] \( *)
(*\*SubsuperscriptBox[\(\[Integral]\), \(\(-\[Pi]\)/q\), \(\[Pi]/q\)]*)
(*\*FractionBox[\(\[DifferentialD]*)
(*\*SubscriptBox[\(k\), \(1\)]\), \(2  \[Pi]/q\)] *)
(*\*SuperscriptBox[\(E\), \(I \(( *)
(*\*SubscriptBox[\(k\), \(1\)] *)
(*\*SubscriptBox[\(x\), \(1\)] + *)
(*\*SubscriptBox[\(k\), \(2\)] *)
(*\*SubscriptBox[\(x\), \(2\)])\)\)] *)
(*\*SubscriptBox[\(\[Psi]\), \(r\)] \(( *)
(*\*SubscriptBox[\(k\), \(1\)], *)
(*\*SubscriptBox[\(k\), \(2\)])\)\)\)\),*)


(* ::Text:: *)
(*with Subscript[\[Psi], r](Subscript[k, 1],k2)\[Congruent] \[Psi](Subscript[k, 1]+2\[Pi] r\[CenterDot]p/q,Subscript[k, 2]).*)
(**)
(*As a result, the iterative Schrodinger E.q. gives the so-called iterative Harper Equation:*)


(* ::DisplayFormulaNumbered:: *)
(*-t(E^(I Subscript[k, 2]) Subscript[\[Psi], r+1](Subscript[k, 1],Subscript[k, 2])+E^(-I Subscript[k, 2]) Subscript[\[Psi], r-1](Subscript[k, 1],Subscript[k, 2])+2 cos(Subscript[k, 1]+2\[Pi] r\[CenterDot]p/q)Subscript[\[Psi], r](Subscript[k, 1],Subscript[k, 2]))=Subscript[E, r] Subscript[\[Psi], r](Subscript[k, 1],Subscript[k, 2]) .*)


(* ::Section:: *)
(*Hofstadter Butterfly*)


ClearAll[HarperMatrix];
(*iterative Schrodinger equation for Hofstadter Model (known as Harper Equation)*)
HarperMatrix[{p_,q_},nu_:0]:=Module[{r},
    N@SparseArray[{{r_,r_}->2Cos[2Pi*r*p/q + nu], {i_,j_}/;Abs[i-j]==1->1},{q,q}]
]
(*HarperMatrix[{2,7}]//MatrixForm*)

ClearAll[pts,data];
pts[N_]:=Block[{eigvals,pqPairs,fractionlist},
	pqPairs=Subsets[Range@N,{2}]~Select~(#[[1]]<=#[[2]]&);(*here {p,q} does not have to be coprime pairs*)
	fractionlist=DeleteDuplicates[Divide@@@pqPairs];
	pqPairs=({Numerator@#, Denominator@#}&/@fractionlist);
	eigvals= Eigenvalues/@(HarperMatrix[#,0]&/@pqPairs);
	
	(Table[{#[[1]],x},{x,#[[2]]}]&/@MapThread[{#1,#2}&,{fractionlist,eigvals}])~Flatten~1
];

(data=pts[30];)//AbsoluteTiming


ListPlot[data,Frame->True,Axes->True,PlotRange->{{0,1},Automatic},FrameLabel->{"Flux/Plaquette","\!\(\*SubscriptBox[\(E\), \(n\)]\)"},FrameTicks->{{Automatic,{}},{{{0,"-2\[Pi]"},{1/4,"-\[Pi]"},{1/2,"0"},{3/4,"\[Pi]"},{1,"2\[Pi]"}},{}}},ImageSize->500]

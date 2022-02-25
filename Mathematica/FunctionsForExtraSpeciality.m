(* ::Package:: *)

(* ::Subsection::Initialization:: *)
(*(*Ini*)*)


(* ::Input::Initialization:: *)
(*
SetDirectory["Z:/Documents/math/Howard"];
<<FunctionsMatrixAnalysis.m;
<<FunctionsComplexAnalysis.m;
<<FunctionsFiniteGroups.m;
<<FunctionsGF2.m;
*)
(* SetDirectory["Z:/Documents/math/Howard/BFstuff"]; *)

dB = 10*Log[10,N[#]]&;
Clear[mm,rr]


(* ::Subsubsection::Initialization::Closed:: *)
(*(*ChordalDistance[C1, C2] - actually squared chordal distance ... *)*)


(* ::Input::Initialization:: *)
(* The assumption is that the codewords are tall codeowrd matrices with Frobenius norm 1 *)



(* ::Input::Initialization:: *)
ChordalDistance = Module[{C1,C2,len,innr},
C1=#1; C2 = #2;
(* Construct projectors from the codewords*)
len = Length[Transpose[C1]];
innr = len*HermConj[C1] . C2;
len- Tr[HermConj[innr] . innr]
]&;


(* ::Input:: *)
(*ChordalDistance = Module[{C1,C2,len,P1,P2},*)
(*C1=#1; C2 = #2;*)
(*(* Construct projectors from the codewords*)*)
(*len = Length[Transpose[C1]];*)
(*P1 = len*C1 . HermConj[C1];P2 = len*C2 . HermConj[C2];*)
(*Tr[(P1-P2) . (P1-P2)]/2*)
(*]&;*)


(* ::Subsubsection::Initialization::Closed:: *)
(*(*PauliPicker[{a, b}]  - picks a Pauli matrix determined by binary vector: Hermitianized Weyl matrix *)*)


(* ::Input::Initialization:: *)
(* to avoid confusion, sigy changed to correspond to Ashikhmin2010. This is just the choice of sign Sqrt[omega] = +-I, which then pops up in commutation relations *)
sigx={{0,1},{1,0}}; sigz={{1,0},{0,-1}}; sigy = -{{0,-I},{I,0}};


(* ::Input:: *)
(*(* Weyl picker is directly defined algebraicly, as powers of sigx and sigz *)*)


(* ::Input::Initialization:: *)
PauliPicker = Module[{inx},
inx=FromDigits[#1,2];If[inx==0,unit2,If[inx==1,sigz,If[inx==2,sigx,sigy]]]
]&;

WeylPicker = Module[{inx1,inx2},
inx1=#1[[1]]; inx2 = #1[[2]]; If[inx1==0,unit2,sigx] . If[inx2==0,unit2,sigz] 
]&;
(* Adding the symplectic inner product of a,b one gets the Hermitian version, which is just the same as Pauli Picker *)
HermWeylPicker = Module[{inx1,inx2},
inx1=#1[[1]]; inx2 = #1[[2]]; I^(-inx1*inx2)*If[inx1==0,unit2,sigx] . If[inx2==0,unit2,sigz] 
]&;

tups = Tuples[{0,1},{2}]
mfl[PauliPicker[#] &/@ tups]
mfl[WeylPicker[#] &/@ tups]
mfl[HermWeylPicker[#] &/@ tups]
(* defining relation of Pauli Picker: (1,1) gives sigy, which is -I sigx sigz *)
first = {1,0}; second =  {0,1}; SymplInnerNotMod2[first,second]
mfl[{I^-SymplInner[first,second]*PauliPicker[first] . PauliPicker[second],PauliPicker[first+second],-I*sigx . sigz,sigy}]



(* ::Input::Initialization:: *)
(* Dual Pauli Picker, takes powers of Y instead of X *)
PauliPickerYZ = Module[{inx},
inx=FromDigits[#1,2];If[inx==0,unit2,If[inx==1,sigz,If[inx==2,sigy,-sigx]]]
]&;

first = {1,0}; second =  {0,1}; SymplInnerNotMod2[first,second]
mfl[{I^-SymplInner[first,second]*PauliPickerYZ[first] . PauliPickerYZ[second],PauliPickerYZ[first+second],-I*sigy . sigz,-sigx}]


(* ::Input:: *)
(*(* the albgebra of the X-matrices, up to sign, is given by I^SymplInnr and the binary sum *)*)
(*bits1 = {0,1}; *)
(*bits2 = {1,0};*)
(*m1 = PauliPicker[bits1];m2 = PauliPicker[bits2];*)
(*md=Mod[bits1 . bits2,2]*)
(*si = SymplInner[bits1,bits2]*)
(*b12 = bits1[[1]]*bits2[[2]]*)
(* {m1 . m2, (-1)^b12*I^si*PauliPicker[Mod[bits1+bits2,2]]} //mfl*)


(* ::Subsubsection::Initialization::Closed:: *)
(*(*TensorPauliPicker[{avec, bvec}] - picks a 2^m dim matrix determined by binary 2m vector*)*)


(* ::Input::Initialization:: *)
(* this handles any 2^m dim *)
TensorPauliPicker = Module[{abvec,len,aa,bb,abpairs,ComponentMatrices},
abvec=#1;
len = Length[abvec]/2;
aa = Take[abvec,len];
bb = Take[abvec,-len];
abpairs = Transpose[{aa,bb}];
ComponentMatrices = PauliPicker[#] &/@ abpairs;
TensorMany[ComponentMatrices]
]&;

TensorPauliPicker[{1,0,0,1,0,1}] //mf
TensorPauliPicker[{1,0,0,1,1,0}] //mf
TensorPauliPicker[{1,0,0,1,0,0}] //mf


(* ::Input::Initialization:: *)
(* this picks from matrices Y and Z instead of X,Z *)
TensorPauliPickerYZ = Module[{abvec,len,aa,bb,abpairs,ComponentMatrices},
abvec=#1;
len = Length[abvec]/2;
aa = Take[abvec,len];
bb = Take[abvec,-len];
abpairs = Transpose[{aa,bb}];
ComponentMatrices = PauliPickerYZ[#] &/@ abpairs;
TensorMany[ComponentMatrices]
]&;

{TensorPauliPicker[{1,0,0,1,0,1}],TensorPauliPickerYZ[{1,0,0,1,0,1}]} //mfl
{TensorPauliPicker[{1,0,0,1,1,0}],TensorPauliPickerYZ[{1,0,0,1,1,0}] } //mfl
{TensorPauliPicker[{1,0,0,1,0,0}] ,TensorPauliPickerYZ[{1,0,0,1,0,0}]} //mfl


(* ::Subsubsection::Initialization:: *)
(*(*MakePauli[avec, bvec] - the same as above, two arguments*)*)


(* ::Input::Initialization:: *)
MakePauli=Module[{avec,bvec},avec = #1; bvec = #2;TensorPauliPicker[Join[avec,bvec]]]&; 


(* ::Subsection::Initialization::Closed:: *)
(*(*InitializeSp[mm] : Initialize structures for Sp (2 mm, 2): WHn, WHnDI, Zerom, Omega, AllBinFullRank, AllBinSymm, HWmodCenter, HWgroup, HWdiag*)*)


(* ::Input::Initialization:: *)
InitializeSp=Module[{},
mm=#1;
Nn = 2^mm;
WHn = WHmatrix[Nn]/Sqrt[Nn]; 
WHnDI = WHmatrixDI[Nn]/Sqrt[Nn];
Zerom = ZeroMatrix[mm]; unitm = IdentityMatrix[mm];  unitN = IdentityMatrix[Nn];
Omega = BlockMatrix[Zerom,unitm,unitm,Zerom];

(* One generate GL(n,2) *)
AllBinFullRank = GenerateGLn2[mm];
          (* Generate additive group of symmetric m x m matrices *)
AllBinSymm =GenerateSn2[mm];

(* Order of Symm(m,2) *)
SymmBinOr = 2^(mm (mm+1)/2);
(* Order of invertible binary matrices *)
InvBinOr = Product[2^mm - 2^ii,{ii,0,mm-1}];
If[SymmBinOr-Length[AllBinSymm]!=0, Print["Error in generating all binary symmetric!?"]];
If[InvBinOr-Length[AllBinFullRank]!=0, Print["Error in generating all binary full rank!?"]];
(* Order of Sp(2m,2) *)
SpOr = 2^(mm^2)  Product[2^(2*ii) - 1,{ii,1,mm}];
(* Grassmannian dimensions *)
BinGrassDims =QBinomial[mm,#,2] &/@ (Range[mm+1]-1);

(* the binary representation of the matrix indexes *)
(* these are not needed: MatInxBin = Table[Reverse[IntegerDigits[jj,2,mm]],{jj,0,Nn-1}]
MatInxBinRev = Reverse[#]&/@ MatInxBin; *)
InxBin = Table[IntegerDigits[jj,2,mm],{jj,0,Nn-1}];
(* Heisenberg-Weyl group modulo center  I^m *)
(* It is important that these are in the order of tups.  ALso, these Tups should not be changed. This is used in inverse operations. *)
TheTupsHW = Tuples[{0,1},2mm]; Length[TheTupsHW];
HWmodCenter = TensorPauliPicker[#] &/@ TheTupsHW; 
HWgroup = (ter = #; Sequence@@(#*ter &/@ {1,I,-1,-I})) &/@ HWmodCenter; 
(* the diagonal part of HW-group *)
HWdiag = Complement[If[Norm[#-DiagonalMatrix[Diagonal[#]]]<0.0001,#] &/@ HWmodCenter,{Null}]; 
HWdiag = (ter = #; Sequence@@(#*ter &/@ {1,I,-1,-I})) &/@ HWdiag; 
]&;


(* ::Subsection::Initialization::Closed:: *)
(*(*InitializeBare[mm] : Initialize WHn, WHnDI, Zerom, Omega,  HWmodCenter, HWgroup, HWdiag*)*)


(* ::Input::Initialization:: *)
InitializeBare=Module[{},
mm=#1;
Nn = 2^mm;
WHn = WHmatrix[Nn]/Sqrt[Nn]; 
WHnDI = WHmatrixDI[Nn]/Sqrt[Nn];
Zerom = ZeroMatrix[mm]; unitm = IdentityMatrix[mm];  unitN = IdentityMatrix[Nn];
Omega = BlockMatrix[Zerom,unitm,unitm,Zerom];

(* the binary representation of the matrix indexes *)
(* these are not needed: MatInxBin = Table[Reverse[IntegerDigits[jj,2,mm]],{jj,0,Nn-1}]
MatInxBinRev = Reverse[#]&/@ MatInxBin; *)
InxBin = Table[IntegerDigits[jj,2,mm],{jj,0,Nn-1}];
(* Heisenberg-Weyl group modulo center  I^m *)
(* It is important that these are in the order of tups.  ALso, these Tups should not be changed. This is used in inverse operations. *)
TheTupsHW = Tuples[{0,1},2mm]; Length[TheTupsHW];
HWmodCenter = TensorPauliPicker[#] &/@ TheTupsHW; 
HWgroup = (ter = #; Sequence@@(#*ter &/@ {1,I,-1,-I})) &/@ HWmodCenter; 
(* the diagonal part of HW-group *)
HWdiag = Complement[If[Norm[#-DiagonalMatrix[Diagonal[#]]]<0.0001,#] &/@ HWmodCenter,{Null}]; 
HWdiag = (ter = #; Sequence@@(#*ter &/@ {1,I,-1,-I})) &/@ HWdiag; 
]&;


(* ::Subsection::Initialization::Closed:: *)
(*(*InitializeVeryBare[mm] : Initialize WHn, WHnDI, Zerom, Omega*)*)


(* ::Input::Initialization:: *)
InitializeVeryBare=Module[{},
mm=#1;
Nn = 2^mm;
WHn = WHmatrix[Nn]/Sqrt[Nn]; 
WHnDI = WHmatrixDI[Nn]/Sqrt[Nn];
Zerom = ZeroMatrix[mm]; unitm = IdentityMatrix[mm];  unitN = IdentityMatrix[Nn];
Omega = BlockMatrix[Zerom,unitm,unitm,Zerom];

(* the binary representation of the matrix indexes *)
(* these are not needed: MatInxBin = Table[Reverse[IntegerDigits[jj,2,mm]],{jj,0,Nn-1}]
MatInxBinRev = Reverse[#]&/@ MatInxBin; *)
InxBin = Table[IntegerDigits[jj,2,mm],{jj,0,Nn-1}];
]&;


(* ::Subsection::Initialization::Closed:: *)
(*(*CreateAllPermMats[m]*)*)


(* ::Input::Initialization:: *)
(* This creates all matrices that permute two tensor product spaces *)


(* ::Input::Initialization:: *)
CreateAllPermMats = Module[{mm,pairs,strt,PermMats},
mm=#1;
pairs = Subsets[Range[mm],{2}];
strt = Table[unit2,{mm}];
PermMats = (ThePair = #;Plus@@((matnow = #; nows = strt; (nows[[#]] = matnow)&/@ ThePair; TensorMany[nows]) &/@ {unit2,sigx,sigy,sigz} )/2) &/@ pairs; 
PermMats
]&;

CreateAllPermMats[3] //mfl


(* ::Subsubsection::Closed:: *)
(*Permutation matrix literally permutes the tensor product spaces *)


(* ::Input:: *)
(*PermMat = Plus@@(Tensor[#,#] &/@ {unit2,sigx,sigy,sigz})/2; PermMat //mf*)
(*Clear[a,b,c,d];Tst  = {{a,b},{c,d}}; Tst2  = {{e,f},{g,h}};*)
(*{on=Tensor[Tst,Tst2],to=PermMat . Tensor[Tst2,Tst] . PermMat,on-to} //mfl*)


(* ::Section:: *)
(*Functions for Clifford group < = > Symplectic mapping*)


(* ::Input::Initialization:: *)
WH2i = WHmatrixDI[2]/Sqrt[2]; (* diagonal prop to identity. This is good if one wants to have the diagonal as close to identity as possible *)
WH2 = WHmatrix[2]/Sqrt[2]; (* diagonal [1,-1]. This is useful if one wants to have symmetric H, leading to Hermitian matrices.  *)



(* ::Subsection::Initialization::Closed:: *)
(*(*UnitaryTransformQuotientRep[Sp]*)*)


(* ::Input:: *)
(*(* this function creates a unitary representative of the quotient which is centered at identity *)*)
(*(* Note that WH-matrix with diagonal close to I used *)*)


(* ::Input::Initialization:: *)
UnitaryTransformQuotientRep = Module[{Pp,Sr,Io,Mm,Ss},
Spmat=#1; 
{Pp,Sr,Io,Mm,Ss} = DecomposeSymplectic[Spmat];
(*{UnitaryF1rep[Pp],UnitaryF2rep[Sr],UnitaryFulrep[Io], UnitaryF1rep[Mm],UnitaryF2rep[Ss]}*)
UnitaryF1rep[Pp] . UnitaryF2rep[Sr] . UnitaryFulrepI[Io] . ConjugateTranspose[UnitaryF1rep[Pp] . UnitaryF2rep[Sr]]
]&;


(* ::Subsection::Initialization::Closed:: *)
(*(*Unitary representations of the Sp subgroups*)*)


(* ::Input:: *)
(*(* The example of F2 type from Calderbank book draft, p. 31 *)*)
(*P={{0,1,0},{1,0,1},{0,1,1}};*)
(*(v=#; I^(v . P . v)) &/@ Table[Reverse[IntegerDigits[jj,2,3]],{jj,0,7}]*)


(* ::Input::Initialization:: *)
(* this version works *)
UnitaryF1rep = Module[{Pp,inxres,RepMat},
Pp=#1;
(* the resulting basis vector from permutation *)
inxres = (FromDigits[#,2]+1) &/@ Mod[InxBin . Pp,2];
RepMat = Table[0,{Nn},{Nn}];
(RepMat[[#,inxres[[#]]]] = 1) &/@ Range[Nn];
RepMat
]&;

(* directly in terms of Z-matrices and units, after decomposing symmetric to set of vector outer products *)
UnitaryF2repZ = Module[{symmi,AllVecs,diag,vece,zedPows},
symmi = #1;
diag = If[Norm[symmi]>0,
AllVecs = DecomposeBinSymmToOuterProducts[symmi];
Times@@((vece  =#;
zedPows = {-1,1}^# &/@ vece;
(*Diagonal of Operator creating the automorphism for this vector *)
Table[1,{Nn}] - I* TensorMany[zedPows]) &/@ AllVecs)/Sqrt[2]^Length[AllVecs],
(* else if the symmteric matrix is zero, output identity *)
Table[1,{Nn}]];
DiagonalMatrix[diag]
] &;

(* This works: binary labeling of dimensions *reversed* in F2 !!! *)
UnitaryF2rep = Module[{Pp,v,diag},
Pp=#1;
(* the resulting sign operation *)
diag = (Mod[# . Pp . #,4]) &/@ InxBin;
DiagonalMatrix[(-I)^diag]
]&;

(* This works: binary labeling of dimensions *reversed* in F3 !!! *)
UnitaryF3rep = Module[{Pp,v,diag},
Pp=#1;
(* the resulting sign operation *)
diag = (Mod[# . Pp . #,4]) &/@ InxBin;
WHn . DiagonalMatrix[(-I)^diag] . WHn
]&;

(* This works: Reversed so that symplectic identity leads to unitary identity *)
UnitaryFulrep = Module[{TheDiag,mats},
TheDiag = Take[Diagonal[#1],mm];
mats = If[#==1,unit2,WH2] &/@ TheDiag;
If[Length[TheDiag]>1,KroneckerProduct[Sequence@@mats],Sequence@@mats]
]&;

(* This works: Reversed so that symplectic identity leads to unitary identity, Use WH-matrix which has diagonal close to identity  *)
UnitaryFulrepI = Module[{TheDiag,mats},
TheDiag = Take[Diagonal[#1],mm];
mats = If[#==1,unit2,WH2i] &/@ TheDiag;
If[Length[TheDiag]>1,KroneckerProduct[Sequence@@mats],Sequence@@mats]
]&;



(* ::Input:: *)
(*(**)
(*UnitaryF1repWrong = Module[{Pp,inxres,RepMat},*)
(*Pp=#1;*)
(*(* the resulting basis vector from permutation *)*)
(*inxres = (FromDigits[Reverse[#],2]+1) &/@ Mod[MatInxBin.Pp,2];*)
(*RepMat = Table[0,{Nn},{Nn}];*)
(*(RepMat[[#,inxres[[#]]]] = 1) &/@ Range[Nn];*)
(*Transpose[RepMat]*)
(*]&;*)
(**)
(*(* this version works *)*)
(*UnitaryF1repOle = Module[{Pp,inxres,RepMat},*)
(*Pp=#1;*)
(*(* the resulting basis vector from permutation *)*)
(*inxres = (FromDigits[#,2]+1) &/@ Mod[MatInxBinRev.Pp,2];*)
(*RepMat = Table[0,{Nn},{Nn}];*)
(*(RepMat[[#,inxres[[#]]]] = 1) &/@ Range[Nn];*)
(*RepMat*)
(*]&;*)
(**)
(*(* This works: binary labeling of dimensions *reversed* in F2 !!! *)*)
(*UnitaryF2repOle = Module[{Pp,v,diag},*)
(*Pp=#1;*)
(*(* the resulting sign operation *)*)
(*diag = (v=Reverse[#]; Mod[v.Pp.v,4]) &/@ MatInxBin;*)
(*DiagonalMatrix[(-I)^diag]*)
(*]&;*)
(**)
(*(* This works: binary labeling of dimensions *reversed* in F3 !!! *)*)
(*UnitaryF3repOle = Module[{Pp,v,diag},*)
(*Pp=#1;*)
(*(* the resulting sign operation *)*)
(*diag = (v=Reverse[#]; Mod[v.Pp.v,4]) &/@ MatInxBin;*)
(*WHn.DiagonalMatrix[(-I)^diag].WHn*)
(*]&;*)
(**)
(*(* This generalizes the ul-matrices in Calderbank et al, any order of tensor products possible. This is a subgroup *)*)
(*UnitaryFulrepWrong = Module[{TheDiag,mats},*)
(*TheDiag = Take[Diagonal[#1],mm];*)
(*mats = If[#\[Equal]1,WH2,unit2] &/@ TheDiag;*)
(*If[Length[TheDiag]>1,KroneckerProduct[Sequence@@mats],Sequence@@mats]*)
(*]&;*)
(**)*)


(* ::Input:: *)
(*UnitaryF1rep[#] &/@ AllBinFullRank//mfl;*)
(*UnitaryF2rep[#] &/@ AllBinSymm//mfl;*)
(*UnitaryF3rep[#] &/@ AllBinSymm//mfl;*)
(*UnitaryFulrep[#] &/@ Fulset //mfl;*)


(* ::Subsubsection::Closed:: *)
(*Check Robert' s direct mapping of F2 - subgroup vs.  Z - mapping*)


(* ::Input:: *)
(*(* Take Random symmetric matrix, make corresponding symplectic  *)*)
(*tst =AllBinSymm[[RandomInteger[{1,Length[AllBinSymm]}]]];*)
(*SpS = SymplecticF2FromSymm[tst];*)
(*(* Create a unitary rep either by using Robert's prescription, or though HW Z-matrices *)*)
(*SpU=UnitaryF2rep[tst];*)
(*SpUz = UnitaryF2repZ[tst];*)
(*(* These are not necessarily similar *)*)
(*prod=HermConj[SpU] . SpUz; Conjugate[prod[[1,1]]]*prod //mf*)
(**)
(*(* Action of both unitaries on full HW group *)*)
(*(* if the output is 0, adjoint action of the unitary rep gives precisely the Sp-action of the corresponding amtrix, if 2, one gets an additional minus-sign *)*)
(*(* all HW indices *)*)
(*tups = Tuples[{0,1},{2mm}];*)
(*(tup = #;*)
(*HWmat = TensorPauliPicker[tup]; *)
(*dir =TensorPauliPicker[Mod[tup . SpS,2]];*)
(*adj=ConjugateTranspose[SpU] . HWmat . SpU; *)
(*adjZ=ConjugateTranspose[SpUz] . HWmat . SpUz; *)
(*{If[Norm[adj-dir]==0 ,"+",If[Norm[adj+dir]==0,"-",0]],If[Norm[adjZ-dir]==0 ,"+",If[Norm[adjZ+dir]==0,"-",0]] } ) &/@ tups*)


(* ::Subsubsection::Closed:: *)
(*Check Adjoint action of all rank - 1  symmetric matrices. All have sign uncertainty. CHECK where this appears from!*)


(* ::Input:: *)
(*mm=3; vecs = Drop[Tuples[{0,1},{mm}],1]*)
(*(* all symmetric matrices that are outerproduct of one vector *)*)
(*SymmMatsFromVecs = Mod[(vec = #; vec  =  {#} &/@ vec; (vec . Transpose[vec])) ,2] &/@ vecs;*)
(**)
(*(* the fully occupied one is not identified as rank-1*)*)
(*tst = SymmMatsFromVecs[[7]]; tst //mf*)
(*(* and accordingly the resulting representation matrix is not of the form 1 + Z !!?? *)*)
(*UnitaryF2repZ[tst] //mf*)
(**)
(*(* The adjoint action of the unitary representatives. All actions ahve the sign uncertainty. *)*)
(*(* if the output is 0, adjoint action of the unitary rep gives precisely the Sp-action of the corresponding amtrix, if 2, one gets an additional minus-sign *)*)
(*(* all HW indices *)*)
(*tups = Drop[Tuples[{0,1},{2mm}],1];*)
(**)
(*(SymMat = #;*)
(*SpS = SymplecticF2FromSymm[SymMat];*)
(*SpUz = UnitaryF2repZ[SymMat];*)
(*Union[(tup = #;*)
(*HWmat = TensorPauliPicker[tup]; *)
(*dir =TensorPauliPicker[Mod[tup . SpS,2]];*)
(*adj=ConjugateTranspose[SpUz] . HWmat . SpUz; *)
(*If[Norm[adj-dir]==0 ,"+",If[Norm[adj+dir]==0,"-",0]] ) &/@ tups]*)
(*) &/@ SymmMatsFromVecs*)


(* ::Subsection::Initialization::Closed:: *)
(*(*UnitaryTransformFromSymplectic[Sp]*)*)


(* ::Input::Initialization:: *)
UnitaryTransformFromSymplectic = Module[{Pp,Sr,Io,Mm,Ss},
Spmat=#1;
{Pp,Sr,Io,Mm,Ss} = DecomposeSymplectic[Spmat];
(*{UnitaryF1rep[Pp],UnitaryF2rep[Sr],UnitaryFulrep[Io], UnitaryF1rep[Mm],UnitaryF2rep[Ss]}*)
UnitaryF1rep[Pp] . UnitaryF2rep[Sr] . UnitaryFulrep[Io] . UnitaryF1rep[Mm] . UnitaryF2rep[Ss]
]&;


(* ::Input:: *)
(*(* The inner automorphism group, with center *)*)
(*InnerAuts= GenerateGroup[Join[{(1-I)/Sqrt[2] HWmodCenter[[1]]},HWgroup],5];*)


(* ::Input:: *)
(*SpS = Mod[tsts[[RandomInteger[len]+1]] . tsts[[RandomInteger[len]+1]] . tsts[[RandomInteger[len]+1]],2]; SpS//mf;*)
(*SpU = UnitaryTransformFromSymplectic[SpS]; Sqrt[8]*SpU //mf;*)
(**)
(*SpS2=Mod[tsts[[RandomInteger[len]+1]] . tsts[[RandomInteger[len]+1]] . tsts[[RandomInteger[len]+1]],2]; *)
(*SpU2 = UnitaryTransformFromSymplectic[SpS2]; Sqrt[8]*SpU2 //mf;*)
(**)
(*(* This is really a group homomorphism, modulo inner automorphisms. *)*)
(*ClassSp = UnitaryTransformFromSymplectic[SpS . SpS2] . # &/@ InnerAuts;*)
(*ClassDir = SpU . SpU2 . # &/@ InnerAuts; Complement[ClassSp,ClassDir] //Length*)


(* ::Input:: *)
(*(* Action on full HW-group, paramterized by the different ab-vectors. The results are indeed + or  - times the directly transformed *)*)
(*Union[(tup = #; adj=ConjugateTranspose[SpU] . TensorPauliPicker[tup] . SpU; dir =TensorPauliPicker[Mod[tup . SpS,2]];*)
(*If[Norm[adj-dir]==0 ,"+",If[Norm[adj+dir]==0,"-",0]]*)
(* ) &/@ Tuples[{0,1},2 mm]]*)


(* ::Subsubsection::Closed:: *)
(*Trying to find the unitary corresp to a Sympl transform directly*)


(* ::Input:: *)
(*SpU=Plus@@Table[Subscript[z, kk]HWmodCenter[[kk]],{kk,1,Nn^2}]; *)
(*SpS = F2set[[3]]; SpS //mf*)
(*(**)
(*Union[(tup = #; dir =TensorPauliPicker[Mod[tup.SpS,2]];ShouldBeId=dir.HermConj[SpU].TensorPauliPicker[tup].SpU; *)
(*Tr[ShouldBeId] //simpconjus//Expand*)
(* ) &/@ Tuples[{0,1},2 mm]]*)
(**)*)


(* ::Subsection::Initialization::Closed:: *)
(*(*UnitaryTransformFromSymplecticZ[Sp]*)*)


(* ::Input:: *)
(*(* Using the 1+i Z form of F2 representation AND the WH-matrix with diagonal proportional to identity *)*)


(* ::Input::Initialization:: *)
UnitaryTransformFromSymplecticZ = Module[{Pp,Sr,Io,Mm,Ss},
Spmat=#1; 
{Pp,Sr,Io,Mm,Ss} = DecomposeSymplectic[Spmat];
(*{UnitaryF1rep[Pp],UnitaryF2rep[Sr],UnitaryFulrep[Io], UnitaryF1rep[Mm],UnitaryF2rep[Ss]}*)
UnitaryF1rep[Pp] . UnitaryF2repZ[Sr] . UnitaryFulrepI[Io] . UnitaryF1rep[Mm] . UnitaryF2repZ[Ss]
]&;


(* ::Input:: *)
(*SpS = Mod[tsts[[RandomInteger[len]+1]] . tsts[[RandomInteger[len]+1]] . tsts[[RandomInteger[len]+1]],2]; Spmat //mf*)
(*SpU = UnitaryTransformFromSymplecticZ[SpS]; Sqrt[8]*SpU //mf*)


(* ::Input:: *)
(*(* Action on full HW-group, paramterized by the different ab-vectors. The results are indeed + or  - times the directly transformed *)*)
(*Union[(tup = #; adj=ConjugateTranspose[SpU] . TensorPauliPicker[tup] . SpU; dir =TensorPauliPicker[Mod[tup . SpS,2]];*)
(*If[Norm[adj-dir]==0 ,"+",If[Norm[adj+dir]==0,"-",0]]*)
(* ) &/@ Tuples[{0,1},2 mm]]*)


(* ::Subsection::Initialization::Closed:: *)
(*(*SymplecticTransformFromUnitary[UnitaryMatrix]*)*)


(* ::Input::Initialization:: *)
(* constructs symplectic presentation of unitary matrix, if the matrix belongs to the Unitary rep of Sp. Otherwise, output "False" *)


(* ::Input::Initialization:: *)
SymplecticTransformFromUnitary = Module[{SpU,mm,TheTups,inxs,SpS,origvec,tup,adj,dir,ActionCheck,trans,vec},
SpU= #;
mm = Log[2,Length[SpU]];
TheTups = Tuples[{0,1},2mm];
If[!UnitaryMatrixQ[SpU],Print["Argument Not Unitary matrix!"],
(* The indexes of the vectors that the symplectic identity matrix is transformed to *)
(* matrix elements found from action on columns in unit *)
inxs = (trans = HermConj[SpU] . TensorPauliPicker[#] . SpU; Position[Abs[Tr[ConjugateTranspose[trans] . #]] &/@ HWmodCenter,Nn][[1,1]]) &/@ IdentityMatrix[2*mm];
SpS = Plus@@((inx = inxs[[#]]; origvec = {#} &/@ (IdentityMatrix[2*mm][[#]]); vec = { TheTups[[inx]]}; origvec . vec ) &/@ Range[2*mm]);
(* test action on full HW-group *)
(* Action on full HW-group, paramterized by the different ab-vectors. The results are indeed + or  - times the directly transformed *)
ActionCheck = (tup = #; adj=ConjugateTranspose[SpU] . TensorPauliPicker[tup] . SpU; dir =TensorPauliPicker[Mod[tup . SpS,2]];
If[Norm[adj-dir]==0 ,"+",If[Norm[adj+dir]==0,"-",0]]
 ) &/@ Tuples[{0,1},2 mm];
ActionCheck = Complement[ActionCheck ,{"+","-"}];
If[ActionCheck=={},SpS,False]]
]&;


(* ::Subsection::Initialization::Closed:: *)
(*(*UnitaryComponentsFromSymplectic[Sp]*)*)


(* ::Input::Initialization:: *)
UnitaryComponentsFromSymplectic = Module[{Pp,Sr,Io,Mm,Ss},
Spmat=#1; 
{Pp,Sr,Io,Mm,Ss} = DecomposeSymplectic[Spmat];
(*{UnitaryF1rep[Pp],UnitaryF2rep[Sr],UnitaryFulrep[Io], UnitaryF1rep[Mm],UnitaryF2rep[Ss]}*)
{UnitaryF1rep[Pp],UnitaryF2rep[Sr],UnitaryFulrep[Io],UnitaryF1rep[Mm],UnitaryF2rep[Ss]}
]&;


(* ::Section:: *)
(*Functions for MUBness*)


(* ::Subsection::Initialization::Closed:: *)
(*(*UBcheck[B1, B2]  - Check if two bases are unbiased*)*)


(* ::Input::Initialization:: *)
UBcheck = Module[{fir,sec,len},fir=#1;sec=#2; len=Length[fir];If[Union[Flatten[Abs[HermConj[fir] . sec]-Table[1,{len},{len}]/Sqrt[len]]]=={0},1,0]]&;

UBcheck[(unit2+sig1)/Sqrt[2],(unit2+sig2)/Sqrt[2]]


(* ::Subsection::Initialization::Closed:: *)
(*(*MUBcheck[{B}]  - Check if set of bases are unbiased*)*)


(* ::Input::Initialization:: *)
MUBcheck = Module[{set,fir},set=#1; 
(fir = #; UBcheck[fir,#] &/@ set) &/@ set
]&;


(* ::Input::Initialization:: *)
MUB4try = {{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}},{{1/2,I/2,I/2,-(1/2)},{I/2,1/2,-(1/2),I/2},{I/2,-(1/2),1/2,I/2},{-(1/2),I/2,I/2,1/2}},{{1/2,-(1/2),-(1/2),1/2},{1/2,1/2,-(1/2),-(1/2)},{1/2,-(1/2),1/2,-(1/2)},{1/2,1/2,1/2,1/2}}};
MUBcheck[MUB4try] //mf

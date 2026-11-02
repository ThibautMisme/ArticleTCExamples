\\ \r AaTechnique.gp
\\ \r AaMakdisiMachinery.gp
\\ \r AaPairing.gp
\\install("PicTorsBasis","GGG")
\\install("PicLiftTors","GGGLL")

InitBasisJ2(J)=
{
	\\return(PicTorsBasis(J,2,0));
	return(pictorsbasis(J,2,0,1));
}

InitBasisJ2FF(JFF)=
{
	my(basisNicolas);
	basisNicolas=InitBasisJ2(JFF);
	if(basisNicolas == 0,
		print("We couldn't produce the whole space J[2](Fq).");
		print("Please try again with the same inputs.");
		);
	return(basisNicolas);
}

/*
InitBasisJ2GAdic(J,BasisJ2FF)=
{
	my(BasisJ2GAdic);
	BasisJ2GAdic=vector(#BasisJ2FF,i,0);
	for(i=1,#BasisJ2GAdic,
		BasisJ2GAdic[i]=LiftTorsion(J,BasisJ2FF[i]);
		print("Lift Basis: ",i,"/",2*getg(J));
		);
	return(BasisJ2GAdic);
}*/

MirrorE(AW,e,V)=
{
	my(L,f);
	L=mattranspose(e)*AW*V;
	f=vector(#V[1,],i,Mod(random(2),2))~;
	while(L*f!=Mod(1,2),
		f=vector(#V[1,],i,Mod(random(2),2))~;);
	f=V*f;
	return(f);
}

InitWeilMatrix(JFF,CCurve,Basis,Weildata)=
{
	return(matrix(#Basis,#Basis,i,j,WeilPairingFF(JFF,CCurve,Basis[i],Basis[j],Weildata)));
}

InitWeilMatrixNicolas(JFF,Basis)=
{
        return(matrix(#Basis,#Basis,i,j,WeilPairingFFNicolas(JFF,Basis[i],Basis[j])));
}

WeilDual(AW,e,f,V)=
{
	my(Le,Lf,Kef,KV,K,Vp);
	Le=mattranspose(e)*AW;
	Lf=mattranspose(f)*AW;
	Le=MatrixForm(Le);
	Lf=MatrixForm(Lf);
	Kef=concatEquationMatrix(Le,Lf);
	KV=matker(V~)~;
	K=concatEquationMatrix(Kef,KV);
	Vp=matker(K);
	return(Vp);
}

SymplecticBasisCoef(J,WeilMatrix)=
\\returns the coefficients in the basis used to compute WeilMatrice of a symplectic basis 
{
	my(g,listeE,listeF,V,e,f);
	\\print("optimiser: minimiser les coefs de e et f en fonction de la base sous-jacente: évite par la suite de gros Sumj");
	g=getg(J);
	listeE=vector(g,i,0);
	listeF=vector(g,i,0);
	V=Mod(matid(2*g),2);
	for(k=1,g,
		e=V[,1];
		listeE[k]=e;
		f=MirrorE(WeilMatrix,e,V);
		listeF[k]=f;
		V=WeilDual(WeilMatrix,e,f,V);
		);
	return(concat(listeE,listeF));
}

export(CoefNonZero)
export(SumJ)
export(AddFlipTM)
export(getW0)
export(getWn)
export(getKn)
export(getd0)
export(getg)
export(DivAddTM)
\\export(DivAdd)
\\export(DivAdd0)
export(DivSubTM)
\\export(DivSub)
export(mulPointWise)
export(getTNicolas)
export(getpeNicolas)
export(getpMak)
export(geteJ)
export(mulpadicElt)

SymplecticBasis(JFF,BasisJ2FF)=		
\\SymplecticBasis(JFF,CCurve,BasisJ2FF,Weildata)=
\\the matrix of the Weil Pairing in that basis is the symplectic matrix [[0,Id],[Id,0]]
{
	my(WeilMatrix,Coef,sympBasis,Pcoef,P);
	\\WeilMatrix=InitWeilMatrix(JFF,CCurve,BasisJ2FF,Weildata);
	WeilMatrix=InitWeilMatrixNicolas(JFF,BasisJ2FF);
	Coef=SymplecticBasisCoef(JFF,WeilMatrix);
	Coef=lift(Coef);
	SympBasis=vector(#BasisJ2FF,i,0);

/*
	for(k=1,#BasisJ2FF,
		Pcoef=Coef[k];
		Pcoef=CoefNonZero(Pcoef);
		P=BasisJ2FF[Pcoef[1]];
		for(i=2,#Pcoef,
			P=SumJ(JFF,P,BasisJ2FF[Pcoef[i]]););
		SympBasis[k]=P;
		);

	A=SympBasis;
*/

	parfor(k=1,#BasisJ2FF,
		Pcoef=Coef[k];
                Pcoef=CoefNonZero(Pcoef);
                P=BasisJ2FF[Pcoef[1]];
                for(i=2,#Pcoef,P=SumJ(JFF,P,BasisJ2FF[Pcoef[i]]););
		P,
			r,SympBasis[k]=r);
	
	\\B=SympBasis;print("A == B ??",A==B);

	return(SympBasis);
}

InitSympBasisJ2FF(JFF)=
{
	my(BasisJ2FF,SBasisJ2FF,MFrobNicolas,MAutsNicolas,MWeilNicolas);
	print("\n");
	print("Beginning Basis J[2](Fq) (Weil Pairing use)");
	[BasisJ2FF,MFrobNicolas,MAutsNicolas,MWeilNicolas]=InitBasisJ2FF(JFF);

	\\print(lift(InitWeilMatrix(JFF,CCurve,BasisJ2FF,Weildata)),"\B");

	\\SBasisJ2FF=SymplecticBasis(JFF,CCurve,BasisJ2FF,Weildata);
	SBasisJ2FF=SymplecticBasis(JFF,BasisJ2FF);
	print("end Basis J[2]", "\n");
	return(SBasisJ2FF);
}



	

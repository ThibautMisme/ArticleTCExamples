\\**************************************************************************************************
\\*                                            COMMUN                                              *
\\**************************************************************************************************



\\********************************************* VRAC ************************************************

IsInList(list,elt)=
{
        my(flag,k);
        flag=0;
        k=1;
        while((flag==0)&&((#list)>=k),
                if(list[k]==elt,flag=1);
                k=k+1;
                );
        return(flag);
}

EquationMatrix(M,dimV,flagF2)=
{
	my(res);
	if(M!=[;],res=matker(M~)~,
		res=matid(dimV);
		if(flagF2,res=Mod(res,2));
		);
	return(res);
}


concatEquationMatrix(K1,K2)=
{
        return(concat(K1~,K2~)~);
}

\\********************************************** DICO TC *********************************************
\\F is for F2^{2g} and is a representation of a TC as vector
\\Nt is for NumberTC and is a represnetation of a TC under the bijection {TC}<->[1,2^{2g}]  ATTENTION: 0<->2^{2G}-1
\\No is for NumberOTC and is a representation of a TC under the bijection {OTC}<->[1,|OTC|]


DicoTCFtoNt(v)=
\\F2^2g -> [0,2^2g-1] ATTENTION: 0<->2^{2G}-1
{
        my(nb);
        nb=0;
        for(i=1,#v,
                nb=nb+lift(v[i])*2^(i-1);
                );
        return(nb);
}

DicoTCNttoF(g,nb)=
\\ATTENTION: 0<->2^{2G}-1 [0,2^2g-1] -> F2^2g
{
        my(v,a);
        a=nb;
        v=vector(2*g,i,Mod(0,2))~;
        for(i=1,2*g,
                v[i]=Mod(a%2,2);
                a=a-lift(v[i]);
                a=a/2;
                );
        return(v);
}

DicoTCNttoNo(nb,listOTCNt)=
\\ATTENTION: 0<->2^{2G}-1 [1,2^2g] -> [1,|OTC|] if nb represents a OTC indeed
{
        my(flag,No);
        flag=0;
        No=0;
        for(i=1,#listOTCNt,
                if(nb==listOTCNt[i],
                        flag=1;
                        No=i;
                        );
                );
        return([flag,No]);
}

DicoTCNotoNt(no,listOTCNt)={return(listOTCNt[no]);}

DicoTCNotoF(no,g,listOTCNt)=
{
	my(Nt,v);
	Nt=DicoTCNotoNt(no,listOTCNt);
	v=DicoTCNttoF(g,Nt);
	return(v);
}

DicoTCFtoNo(v,listOTCNt)=
{
	my(Nt,No);
	Nt=DicoTCFtoNt(v);
	No=DicoTCNttoNo(Nt,listOTCNt);
	if(No[1]==1,No=No[2],print("error: Nt is not an OTC"));
	return(No);
}

\\*************************************** QTF ********************************************

QuadCanonique(g,v)=
\\g = genus = #(ei) = #fi
{       
        my(res);
        res=0;
        for(k=1,g,res=res + v[k]*v[k+g]);
        return(res);
}

QuadB0(B0data,v)=
\\WeilMatrix is supposed to be the symplectic Matrix
{
        \\my(g,WeilMatrix,vB0,ArfB0);
        \\[g,WeilMatrix,vB0,ArfB0]=B0data;
        return(QuadCanonique(B0data[1],v) + mattranspose(B0data[3])*B0data[2]*v);
}



\\********************************************* ***********************************************

PrintOrbitDecomposition(matfact,flagShort)=
{
	my(nbRat,decomp,deg,pow,flag);
	nbRat=0;
	decomp=[];
        for(k=1,#matfact[,1],
                deg=poldegree(matfact[k,1]);
                pow=matfact[k,2];
                flag=1;
                for(l=1,#decomp,
                        if(deg==decomp[l][1],
                                decomp[l][2]=decomp[l][2]+pow;
                                flag=0;
                                );
                        );
                if(flag,
                        decomp=concat(decomp,[[deg,pow]]);
                        );
                );
	if(flagShort==0,print1("Orbit Decomposition: "));for(l=1,#decomp,print1(decomp[l][1],"_",decomp[l][2]," "));print1("\n");
	if(decomp[1][1]==1,
		if(flagShort==0,print("RATIONAL ELEMENT: ",decomp[1][2]));
		nbRat=decomp[1][2];
		);
        return(nbRat);
}


NbRationalCore(Chi,p)=
{
	my(matfact,nbrat);
	if(p==0,matfact=factor(Chi),
		matfact=factorpadic(Chi,p,500));
	nbrat=0;
	for(k=1,#matfact[,1],
		if(poldegree(matfact[k,1])==1,
			nbrat=nbrat+matfact[k,2]));
	return(nbrat);
}

NbRational(ChiOTC,flagTC,p,flagShort)=
{
	my(nbrat);
	nbrat=0;
	if(flagTC,nbrat=NbRationalCore(ChiOTC[1],p)+NbRationalCore(ChiOTC[2],p),
		nbrat=NbRationalCore(ChiOTC,p));
	if(nbrat>=1,
		if(flagShort==0,print("		RATIONAL TC"));
		);
	return(nbrat);
}

NbRationalJ2(ChiJ2,p,flagShort)=
{
        my(nbrat);
        nbrat=NbRationalCore(ChiJ2,p)+1;
        if(nbrat>=2,
                if(flagShort==0,print("         RATIONAL 2-TORSION"));
                );
        return(nbrat);
}

PrintOrbitTC(flagTC,ChiOTC,p,flagShort)=
{
	my(a);
	if(flagTC==0,
		if(p==0,a=PrintOrbitDecomposition(factor(ChiOTC),flagShort),
			a=PrintOrbitDecomposition(factorpadic(ChiOTC,p,500),flagShort)););
	if(flagTC,
		if(p==0,
			a=PrintOrbitDecomposition(factor(ChiOTC[1]),flagShort);
			a=PrintOrbitDecomposition(factor(ChiOTC[2]),flagShort););
		if(p!=0,
			a=PrintOrbitDecomposition(factorpadic(ChiOTC[1],p,500),flagShort);
			a=PrintOrbitDecomposition(factorpadic(ChiOTC[2],p,500),flagShort););
		);
	return(a);
}









\\**************************************************** STEINER COMPLEX ***********************************************

SteinerComplexPoint(g,v,listOTCNt)=
{
	my(Steiner,Nt,votc,Ntsum,couple);
	Steiner=[];
	for(k=1,#listOTCNt,
		Nt=listOTCNt[k];
		votc=DicoTCNttoF(g,Nt);
		Ntsum=DicoTCFtoNt(v+votc);
		couple=Set([Nt,Ntsum]);
		if(IsInList(listOTCNt,Ntsum),
			if(IsInList(Steiner,couple)==0,
				Steiner=concat(Steiner,[couple])));
		);
	return(Steiner);
}


InitSyzygeticLoop(g,listOTCNt,SteinerComplex,nbloop)=
{
	my(Syzy,SC,q0,quadruple);
	\\Syzy=[Set(concat(SteinerComplex[1][1],SteinerComplex[1][2]))];
	Syzy=[];
	for(k=1,nbloop,
		print1(k,", ");
		SC=SteinerComplex[k];
		for(m=1,#SC,
			q0=SC[m];
			for(l=m+1,#SC,
				quadruple=Set(concat(q0,SC[l]));
				if(IsInList(Syzy,quadruple)==0,
					Syzy=concat(Syzy,[quadruple]));
				);
			);
		);
	return(Syzy);
}

InitSteinerComplex(g,listOTCNt)=
{
	my(listeSteiner,v);
	listeSteiner=vector(2^(2*g)-1,i,0);
	for(k=1,#listeSteiner,
		v=DicoTCNttoF(g,k);
		listeSteiner[k]=SteinerComplexPoint(g,v,listOTCNt);
		);
	return(listeSteiner);
}

/*
FromSyzytov(listOTCNt,quadruple)=
{
	my(v,Nt,No);
	v=vector(#listOTCNt,i,Mod(0,2));
	for(k=1,#quadruple,
		Nt=quadruple[k];
		No=DicoTCNttoNo(Nt,listOTCNt);
		v[No]=Mod(1,2);
		);
	return(v);
}
*/

InitSyzygetic(g,listOTCNt,SteinerComplex)=
{
	my(Syzygetic);
	Syzygetic=InitSyzygeticLoop(g,listOTCNt,SteinerComplex,#SteinerComplex);
	return(Syzygetic);
}

FindGeneratorsSyzy(g,listOTCNt,SteinerComplex)=
{
        my(SyzyGen);
        SyzyGen=InitSyzygeticLoop(g,listOTCNt,SteinerComplex,1);
        return(SyzyGen);
}

PrintSyzyGenMagma(SyzyGen)=
{
	print("SyzyGen:=",SyzyGen,";");
	return();
}

\\****************************************************************************************************************

CheckEgalitePointCoordSmooth(PI,QI)=
{
        my(P,Q,flag,lambda);
	if(#PI==2,P=[PI[1],PI[2],1],P=PI);
	if(#QI==2,Q=[QI[1],QI[2],1],Q=QI);
        flag=1;
        lambda=0;
        if(flag,
                for(k=1,3,
                        if(P[k]==0,
                                if(Q[k]!=0,
                                        flag=0;
                                        );
                                );
                        if(P[k]!=0,
                                if(Q[k]==0,
                                        flag=0;
                                        );
                                if(Q[k]!=0,
                                        if(lambda==0,
                                                lambda=P[k]/Q[k];
                                                );
                                        if(lambda!=0,
                                                if((P[k]/Q[k])!=lambda,
                                                        flag=0;
                                                        );
                                                );
                                        );
                                );
                        );
                );
        return(flag);
}

CheckEgalitePointCoord(P,Q)=
{
        my(flag);
	if(default(debug)>1,print("Should use CheckEgalitePointNicolas instead"));
        if((#P==2)||(#P==3),if((#Q==2)||(#Q==3),flag=CheckEgalitePointCoordSmooth(P,Q),flag=0),
                if((#Q==2)||(#Q==3),flag=0,flag=(P==Q)));
        return(flag);
}


install("PtEq","lGGG")
CheckEgalitePointNicolas(CCurve,P,Q)=
{
	return(PtEq(CCurve,P,Q));
}

\\*****************************************************************************************************************
CestNormal()={}

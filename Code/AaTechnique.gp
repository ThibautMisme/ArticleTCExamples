install("PicNorm","GGGGL") 
\\install("JacRed","GL") \\picsetprec
\\install("PicRand","GG") \\min
install("PicRand_THIBAUT","GG")
\\install("PicLift","GGLGL") \\min
\\install("PicTorsPairingInit","GG") \\min
\\install("PicTorsPairing","GGGG")\\min
\\install("FnDiv","GG")
\\install("RiemannRoch","GG")
install("PicEval","GG")
install("ZpXQ_div","GGGGGL")
install("FqM_FqC_mul","GGGG")
install("Fq_mul","GGGG")
install("matkerpadic","GGGGL")
install("FqM_mul","GGGG")
install("ZpXQ_inv","GGGL")

install("Jgetd0","lG")
install("Jgetf","G")
install("Jgetg","lG")
install("Jgetp","G")
install("JgetZ","G")
install("JgetE","G")
install("Jgetpe","G")
install("JgetT","G")
install("JgetW0","G")
install("JgetV","GU")
install("JgetL","G")
install("JgetKV","GU")
install("JgetV","GU")
install("JgetEvalData","G")

\\TECHNIQUE

\\print("WEIL PAIRING: PAS BESOIN DE CALCULER DELTA1 PUISQU'ON PREND QUOTIENT");

getPolCurve(J)={return(Jgetf(J));} \\J[1]
getg(J)={return(Jgetg(J));} \\J[2]
getpMak(J)={return(Jgetp(J));} \\J[6]
getd0(J)={return(Jgetd0(J));} \\J[3]
getdw(J)={return(Jgetd0(J) - Jgetg(J) +1);}
getZ(J)={return(JgetZ(J));} \\J[15]
geteJ(J)={return(JgetE(J));} \\J[7]
getpeNicolas(J)={return(Jgetpe(J));} \\J[8]
getTNicolas(J)={return(JgetT(J));} \\J[5]
getW0(J)={return(JgetW0(J));}\\ J[13]. est le POINT origine
getW1(J)={return(JgetV(J,1));}\\ J[11][1]. est la matrice rpz le point H°(D0)
getW2(J)={return(JgetV(J,2));}\\ J[11][2]
getW1Pol(J)={return(JgetL(J)[1]);}\\ J[4][1]
getW2Pol(J)={return(JgetL(J)[2]);}\\ J[4][2]
getK3(J)={return(JgetKV(J,3));} \\J[12][3]
getWn(J,n)={return(JgetV(J,n));} \\J[11][n]
getKn(J,n)={return(JgetKV(J,n));} \\J[12][n]
getWE(J)={return(JgetEvalData(J));} \\J[14]



\\************************************************* PRINT ************************************************

PrintOrbitDecompositionChi(Chi)=
{
        my(matfact);
        matfact=factor(Chi);
        PrintOrbitDecomposition(matfact);
}

PrintInput(flagTC,flagmethod,flagTFEonephase,flagTH,flagEffective,flagmethodChiJ2,flagmethodTC)=
{
	print("INPUTS:");
	print("		flagTC = ",flagTC);
	print("		flagmethod = ",flagmethod);
	if(flagTC,print("	method TC = ",flagmethodTC));
	print("		flagTFEonephase = ",flagTFEonephase);
	print("		flagTH = ",flagTH);
	print("		flagEffective = ",flagEffective);
	print("		flagmethodChiJ2 = ",flagmethodChiJ2,"\n");
	return();
}


\\********************************************************************************************************




















CurveDiscAffine(PolCurve)=
{
	my(f,v,w,fv,fw,Dv,Dw,Disc);
	f=PolCurve;
	[v,w]=variables(f);
	fv=deriv(f,v);
	fw=deriv(f,w);
	Dv=polresultant(f,fv,v);\\print("Dv = ",Dv);
	Dw=polresultant(f,fw,v);\\print("Dw = ",Dw);
	Disc=polresultant(Dv,Dw,w);

	if(Disc==0,
		Dv=polresultant(f,fv,w);
        	Dw=polresultant(f,fw,w);
        	Disc=polresultant(Dv,Dw,v);
		);
	\\if(Disc==0,
	\\	print("\n","WARNING: CurveDisc HAS FAILED: Dv= ",Dv," Dw = ",Dw,"\n"););

	Disc=abs(Disc);
	return(Disc);
}

BadPrimeDiscAttempt(PolCurve)=
{
	my(fhom,Discz,Discy,Discx,flag,Badprime);
	fhom=polhomogenise(PolCurve,'z);
	Discz=CurveDiscAffine(subst(fhom,'z,1));
	Discy=CurveDiscAffine(subst(fhom,'y,1));
	Discx=CurveDiscAffine(subst(fhom,'x,1));
	flag=1;
	Badprime=[];
	
	if(Discz==0,flag=0);
	if(Discy==0,flag=0);
	if(Discx==0,flag=0);
	if(flag,
		Badprime=factor(Discz)[,1]~;
		Badprime=Set(concat(Badprime,factor(Discy)[,1]~));
                Badprime=Set(concat(Badprime,factor(Discx)[,1]~));
                );
	return([flag,Badprime]);
}

RandomChangeVariable(PolCurve)=
{
	my(coin,u,index,coin2,v,f,coin3);
	coin=random(3)+1;
	u=['x,'y,'z][coin];
	if(coin==1,index=[2,3],
		if(coin==2,index=[1,3],
			index=[1,2]);
		);
	coin2=random(2)+1;
	v=['x,'y,'z][index[coin2]];
	coin3=random(2);
	f=subst(polhomogenise(PolCurve,'z),u,u+((-1)^coin3)*v);
	f=subst(f,z,1);
	return(f);
}

BadPrimeDiscLoop(PolCurve)=
{
	my(PC,attempt,attemptmax,flag,Badprime);
	PC=PolCurve;
	attempt=0;
	attemptmax=30;
	flag=0;
	while((flag==0)&&(attempt<attemptmax),
		attempt=attempt+1;
		[flag,Badprime]=BadPrimeDiscAttempt(PC);\\print("flag = ",flag);
		\\change of variables
		PC=RandomChangeVariable(PC);
		);
	\\if(flag==0,print("ERROR Disc"));
	return([flag,Badprime]);
}

BadPrimeDisc(PolCurve)=
{
	my(flag,Badprime,flagbis,Badprimebis,p,index);
	if(default(debug)>1,print("beginning Badprime computation"));
	[flag,Badprime]=BadPrimeDiscLoop(PolCurve);
	for(k=1,4,
		[flagbis,Badprimebis]=BadPrimeDiscLoop(PolCurve);
		if(flag==0,
			flag=flagbis;
			Badprime=Badprimebis;);
		if(flag,
			if(flagbis,
				index=0;
				for(k=1,#Badprime,
					index=index+1;
					p=Badprime[index];
					if(IsInList(Badprimebis,p)==0,
						if(index==1,Badprime=Badprime[2..#Badprime]);
						if(index==#Badprime,Badprime=Badprime[1..(#Badprime-1)]);
						if((index>1)&&(index<#Badprime),Badprime=concat(Badprime[1..(index-1)],Badprime[(index+1)..#Badprime]));
						index=index-1;
						););););
		);
	if(default(debug)>1,print("end"));
	if(flag==0,print("ERROR Disc"));
	return([flag,Badprime]);
}

/*
BadPrimeDisc(PolCurve)=
{
	return(factor(CurveDisc(PolCurve))[,1]~);
}
*/
/*
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
*/

IsInListNested(list,elt)=
{
	my(flag,k);
	flag=0;
	k=1;
	while((flag==0)&&((#list)>=k),
		if(IsInList(list[k],elt),
			flag=1;);
		k=k+1;
		);
	return(flag);
}

CheckCurveSingular(CCurve,g)=
\\check if CCurve, provided by CrvInit, is smooth or not
\\returns 1 if SINGULAR
{
	\\return((length(CCurve[5])>1));
	return(crvgenus(CCurve)<g); 
}

CheckSingular(PolCurve,g,p)=
{
	my(CC,res);
	CC=crvinit(PolCurve,'t,'a,Mod(1,p));
	if(CC==0,res=1,res=CheckCurveSingular(CC,g));
	return(res);
}

WithdrawSingular(PolCurve,g,Bad)=
{
        my(res,index);
        res=Bad;
        index=0;
        for(k=1,#Bad,
                index=index+1;
                if(CheckSingular(PolCurve,g,Bad[k])==0,
                        print("good p : ",Bad[k]);
                        res=concat(res[1..(index-1)],res[(index+1)..#res]);
                        index=index-1;
                        );
                );
        return(res);
}

CurveBadPrime(PolCurve,g)=
{
	my(flag,Badprime);
	[flag,Badprime]=BadPrimeDisc(PolCurve);print("Bad Primes from Disc: ",Badprime);
        if(flag,
		Badprime=WithdrawSingular(PolCurve,g,Badprime);print("Badprime: ",Badprime);
		);
	return(Badprime);
}

/*
GoodRed(PolCurve,pInitial,BadPrime)=
\\returns p and CCurveFF provided by CrvInit such that the curve has good reduction in p
{
	my(p,flag,CCurveFF);
	p=pInitial;
	if(p==2,
                p=3;
                print("WARNING: pMak=2 not possible. pMak=3 increases");
                );

	\\flag=1;
        \\while(flag,
	while(IsInList(BadPrime,p),
               	p=nextprime(p+1);
               	print("WARNING: pMak is a prime of bad reduction. pMak increases: pMak = ",p);
               	);
		/*
		\\CrvInit badly behave if wild ramification
		CCurveFF=crvinit(Mod(PolCurve,p),'t,'a);
		while(CCurveFF==0,
			p=nextprime(p+1);
                	CCurveFF=crvinit(Mod(PolCurve,p),'t,'a);
			print("CrvInit: Wild Ramification occured");
			);
		
		flag=CheckCurveSingular(CCurveFF);
		if(flag,
			p=nextprime(p+1);
                	print("WARNING: pMak is a prime of bad reduction. pMak increases: pMak = ",p);
			);
		);
	
	CCurveFF=crvinit(Mod(PolCurve,p),'t,'a);
	return([p,CCurveFF]);
}
*/

divpadic(J,a,b)=
{
	return(ZpXQ_div(a,b,getTNicolas(J),getpeNicolas(J),getpMak(J),geteJ(J)));
}

rpzJtoQq2(J,elt)=
\\return the natural representation of elt in Qq given under Nicolas' formalism (in J for instance)
{
	my(a,T,t,pe);
	a=0;
	T=getTNicolas(J);
	t=variable(T);
	pe=getpeNicolas(J);
	for(k=0,poldegree(T),
		a=a + (polcoef(elt,k) +O(pe))*t^k;
		);
	a=Mod(a,T);
	return(a);
	\\return(Mod(elt,getTNicolas(J))*(1+O(getpMak(J)^geteJ(J))));
}

rpzJtoQq(J,elt)=
\\return the natural representation of elt in Qq given under Nicolas' formalism (in J for instance)
{
	return(Mod(elt,getTNicolas(J))*(1+O(getpMak(J)^geteJ(J))));
}

rpzQqtoJ(elt)=
\\get back into the J's formalism from an element provided by rpzJtoQq
{
	my(pol,t,res);
	pol=lift(elt);
	t=variable(pol);
	res=0;
	for(k=0,poldegree(pol),
		res=res + lift(polcoef(pol,k))*t^k;
		);
	return(res);
}
	
rpzJtoQqMatrix(J,M)=
{
	return(matrix(#M[,1],#M[1,],i,j,rpzJtoQq(J,M[i,j])));
}

rpzQqtoJMatrix(M)=
{
	return(matrix(#M[,1],#M[1,],i,j,rpzQqtoJ(M[i,j])));
}


rpzJtoMod(J,elt)=
\\return Mod(elt,T) for elt given under Nicolas' formalism (in J for instance)
{
        return(Mod(elt,getTNicolas(J)));
}

rpzModtoJ(elt)=
\\get back into the J's formalism from an element provided by rpzJtoQq
{
	return(lift(elt));
}

rpzJtoModMatrix(J,M)=
{
        return(matrix(#M[,1],#M[1,],i,j,rpzJtoMod(J,M[i,j])));
}

rpzModtoJMatrix(M)=
{
        return(matrix(#M[,1],#M[1,],i,j,rpzModtoJ(M[i,j])));
}


polcoefMatrix(M,n)=
{
	return(matrix(#M[,1],#M[1,],i,j,polcoef(M[i,j],n)));
}

/*
concatEquationMatrix(K1,K2)=
{
	return(concat(K1~,K2~)~);
}
*/

CheckWinK(J,W,K)=
\\check that W is in KerK
{
	my(M,Mnul);
	\\M=rpzQqtoJMatrix(rpzJtoQqMatrix(J,K)*rpzJtoQqMatrix(J,W));
	M=mulpadicMatrix(J,K,W);
	Mnul=matrix(#M[,1],#M[1,],i,j,0);
	return(M==Mnul);
}


reduceFF(JFF,elt)=
{
        return(Mod(Mod(elt,getpMak(JFF)),Mod(getTNicolas(JFF),getpMak(JFF))));
}
reduceFFCol(JFF,c)=
{
        return(vector(#c,i,reduceFF(JFF,c[i]))~);
}

reduceFFMatrix(JFF,M)=
{
        return(matrix(#M[,1],#M[1,],i,j,reduceFF(JFF,M[i,j])));
}

liftFF(elt)=
{
        return(lift(lift(elt)));
}

liftFFCol(c)=
{
        return(vector(#c,i,liftFF(c[i]))~);
}

liftFFMatrix(M)=
{
        return(matrix(#M[,1],#M[1,],i,j,liftFF(M[i,j])));
}

MulFF(JFF,a,b)=
{
        return(liftFF(reduceFF(JFF,a)*reduceFF(JFF,b)));
}


CheckWinKFF(JFF,W,K)=
{
	my(M,Mnul);
	M=liftFFMatrix(reduceFFMatrix(JFF,K)*reduceFFMatrix(JFF,W));
	Mnul=matrix(#M[,1],#M[1,],i,j,0);
        return(M==Mnul);
}

MulPointWiseFF(JFF,c,W)=
{
	return(matrix(#W[,1],#W[1,],i,j, liftFF(reduceFF(JFF,c[i])*reduceFF(JFF,W[i,j])) ));
}

mulPointWise(J,c,W)=
\\return the pointwise multiplication of the column c with the matrix W
{
	\\return(matrix(#W[,1],#W[1,],i,j, rpzQqtoJ(rpzJtoQq(J,c[i])*rpzJtoQq(J,W[i,j])) ));
	return(matrix(#W[,1],#W[1,],i,j, mulpadicElt(J,c[i],W[i,j])));
}

ColPointWiseFF(JFF,c1,c2)={return(vector(#c1,i,liftFF(reduceFF(JFF,c1[i])*reduceFF(JFF,c2[i])))~);}

/*
RandomPointJFF(JFF)=
{
	my();
	g=getgJ(JFF);
}
*/



RedPointJ(J,W,e)=
\\reduce a point W of J to the precision e
{
        return(matrix(#W[,1],#W[1,],i,j,lift(W[i,j]*(1+O(getpMak(J)^e)))));
}

RedJ(J,e)=
\\reduce the J data to the precision e
{
        \\return(JacRed(J,e));
	return(picsetprec(J,e));
}


InitMakdisiFF(J)=
\\return J for e=1, which is JFF=J(Fq)
{
        return(RedJ(J,1));
}

RandomPointJ(J)=
{
        return(picrand(J,random(50)));
}

RandomPointJFp(J)=
{
	print("long?? Try to INCREASE pInit: It may be possible that J(Fp) is too small to be spanned by the points we have suceeded to compute");
	return(PicRand_THIBAUT(J,random(50)));
}

LiftRandom(J,e,W)=
{
	return(piclift(J,W,e,0,0));
}

InitPairingSettingNicolas(J)=
{
        return(pictorspairinginit(J,2));
}


InitlisteProofFunctionsJ2(J,BasisJ2)=
{
        my(listeProofFunctionsJ2);
        listeProofFunctionsJ2=vector(#BasisJ2,i,ProofFunctionTors(J,BasisJ2[i]));
        return(listeProofFunctionsJ2);
}






DLint(n,g)=
\\return the biggest l s.t. g^l/n
{
	my(x,m);
	x=0;
	m=n;
	while((m%g)!=1,
		x=x+1;
		m=m/g;
		);
	return(x);
}


minv(listev)=
{
	my(gardefou,incrementation,mu,i0,v);
	gardefou=#listev;
	incrementation=1;
	mu=listev[incrementation];
	i0=incrementation;
	while((incrementation<gardefou)&&(mu<0),
		incrementation=incrementation+1;
		mu=listev[incrementation];
        	i0=incrementation;
		);
	if(incrementation==gardefou,
		print("ERROR: D0 Case=1"););
	for(k=incrementation+1,#listev,
		v=listev[k];
		if((v>=0)&&(v<mu),
			mu=v;
			i0=k;
			);
		);
	return([mu,i0]);
}

MatrixForm(V)=
{
	return(matrix(1,#V,i,j,V[j]));
}
ColToMat(listeCol)=
{
	return(matrix(#listeCol[1],#listeCol,i,j,listeCol[j][i]));
}

CoefNonZero(liste)=
{
	my(l);
	l=[];
	for(k=1,#liste,
		if(liste[k]!=0,
			l=concat(l,[k]);););
	return(l);
}

CoefInversible(liste,p)=
{
        my(l);
        l=[];
        for(k=1,#liste,
                if(Mod(liste[k],p)!=0,
                        l=concat(l,[k]);););
        return(l);
}

ChooseColumn(M,s)=
{
	my(liste,dim,index);
	liste=vector(s,i,0);
	dim=#M[1,];
	for(i=1,s,
		index=random(dim)+1;
		while(IsInList(liste,index),
			index=random(dim)+1;);
		liste[i]=index;
		);
	return(liste);
}
			
MatrixFromCol(listeIndex,V)=
{
	return(matrix(#V[,1],#listeIndex,i,j,V[i,listeIndex[j]]));
}

FindSupplementFF(JFF,W,V)=
{
	my(dim,dimV,d,listeIndex,Mp);
	dim=#W[1,];
	dimV=#V[1,];
	d=dimV-dim;
	while(dim!=dimV,
		listeIndex=ChooseColumn(V,d);
		Mp=MatrixFromCol(listeIndex,V);
		dim=matrank(reduceFFMatrix(JFF,concat(W,Mp)));
		);
	return(Mp);
}

PointWise(J,c1,c2)={
	\\return(rpzQqtoJ(rpzJtoQq(J,c1)*rpzJtoQq(J,c2)));
	return(vector(#c1,i, mulpadicElt(J,c1[i],c2[i]) )~);
}

MatrixPointWiseV5FF(JFF,listeIndexV2,listeIndexV3,V2,V3)={return(matrix(#V2[,1],#listeIndexV2,i,j, MulFF(JFF, V2[i,listeIndexV2[j]] , V3[i,listeIndexV3[j]] )));}

FindSupplementV5FF(JFF,W)=
{
	my(V2,V3,dim,dimV5,d,listeIndexV2,listeIndexV3,Mp);
	V2=getWn(JFF,2);
	V3=getWn(JFF,3);
	dim=#W[1,];
	dimV5=5*getd0(JFF) -getg(JFF) +1;
	\\d=dimV5-dim;
	d=getd0(JFF);
	while(dim!=dimV5,
		listeIndexV2=ChooseColumn(V2,d);
		listeIndexV3=ChooseColumn(V3,d);
		Mp=MatrixPointWiseV5FF(JFF,listeIndexV2,listeIndexV3,V2,V3);
		dim=matrank(reduceFFMatrix(JFF,concat(W,Mp)));
		);
	return(Mp);
}

SumVect(c1,c2)={return(vector(#c1,i,c1[i]+c2[i]));}

MulVectFF(JFF,c1,c2)={return(vector(#c1,i,MulFF(JFF,c1[i],c2[i])))~;}

RandomPolCol(J,n)=
{
        my(V,xc,c);
	if(n==0,
		V=getW0(J););
	if(n>0,
		V=getWn(J,n););
        xc=vector(#V[1,],i,random(100))~;
	while(xc==vector(#V[1,],i,0)~,
		xc=vector(#V[1,],i,random(100))~;
		);
        c=V*xc;
        return(c);
}

RandomCol(W)=
{
	my(xc,c);
	xc=vector(#W[1,],i,random(100))~;
	while(xc==vector(#W[1,],i,0)~,
		xc=vector(#W[1,],i,random(100))~;);
        c=W*xc; 
        return(c);
}

RandomPolH(H)=
{
	my(xc,f);
        xc=vector(#H,i,random(20))~;
        while(xc==vector(#H,i,0)~,
                xc=vector(#H,i,random(20))~;);
	f=1;
	for(k=1,#H,
		f=f+xc[k]*H[k];
		);
	return(f);
}

PoltoMakdisiColFF(JFF,poly)=
{
        my(Z,W,Var,gnum,gden);
        Z=getZ(JFF);
        Var=[x,y];
	\\W=vector(#Z,i,liftFF(reduceFF(JFF,substvec(poly,Var,Z[i]))))~;
	W=vector(#Z,i,0);
	gnum=numerator(poly);
        gden=denominator(poly);
        for(i=1,#Z,
        	W[i]=liftFF(reduceFF(JFF,substvec(gnum,Var,Z[i]))/reduceFF(JFF,substvec(gden,Var,Z[i])));
                );
        return(W);
}

rpzJtoJ(J,elt)={return(rpzQqtoJ(rpzJtoQq(J,elt)));}


PoltoMakdisiCore(g,Var,Z)=
{
	my(ll,gy,gynum,gyden,denom,Zg,c);
	ll=matrix(#Z,poldegree(g)+1,i,j,0);
	for(j=0,poldegree(g),
		gy=polcoef(g,j);
		gynum=numerator(gy)*x^j;
		gyden=denominator(gy)*x^0;
		for(i=1,#Z,
                        \\W[i,j]=rpzQqtoJ(rpzJtoQq(J,substvec(gnum,Var,Z[i]))/rpzJtoQq(J,substvec(gden,Var,Z[i])));
                        denom=substvec(gyden,Var,Z[i]);
                        if(Mod(denom,getpMak(J))==0,print("PoltoMakdisi error: division by 0"),
                                ll[i,j+1]=divpadicElt(J,substvec(gynum,Var,Z[i]),denom));
                        );
		);
	Zg=vector(#Z,i,0)~;
	for(i=1,#Z,
		c=0;
		for(j=1,poldegree(g)+1,c=c+ll[i,j]);
		Zg[i]=c;
		);
	return(Zg);
}



PoltoMakdisi(J,Wpol)=
\\Wpol liste de polynomes
{
	my(Z,W,Var,g,gnum,gden,denom,T,pe);
	Z=getZ(J);
	T=getTNicolas(J);
	pe=getpeNicolas(J);
	Z=Mod(Z,T);
	Z=Mod(Z,pe);
	Var=variables(Wpol);
	if(Var==['x],
		Var=concat(Var,['y])
		);
	\\W=matrix(#Z,#Wpol,i,j,rpzQqtoJ(rpzJtoQq(J,substvec(Wpol[j],Var,Z[i]))));
	W=matrix(#Z,#Wpol,i,j,0);
	for(j=1,#Wpol,
		g=Wpol[j];
		g=lift(g);
		gnum=numerator(g);
		Zgnum=PoltoMakdisiCore(gnum,Var,Z);
		gden=denominator(g);
		if(gden==1,W[,j]=Zgnum,
			Zgden=PoltoMakdisiCore(gden,Var,Z);
			for(i=1,#Z,
				denom=substvec(Zgden[i],Var,Z[i]);
				if(Mod(denom,getpMak(J))==0,print("PoltoMakdisi error: division by 0"),
					W[i,j]=divpadicElt(J,substvec(Zgnum[i],Var,Z[i]),denom));
				);
			);
		);
	return(W);
}

/*
CheckEgalitePointCoordSmooth(P,Q)=
{
	my(flag,lambda);
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
	if((#P==2)||(#P==3),if((#Q==2)||(#Q==3),flag=CheckEgalitePointCoordSmooth(P,Q),flag=0),
		if((#Q==2)||(#Q==3),flag=0,flag=(P==Q)));
	return(flag);
}
*/

EvalDivisorNicolas(J,W,f)={return(PicNorm(J,f,vector(#W[,1],i,1)~,W,3));}

ColtoPolFFW2(JFF,c)=
{
        my(V,vc,Vpol,fc);
        V=getWn(JFF,2);
        vc=matsolve(reduceFFMatrix(JFF,V),reduceFFCol(JFF,c));
        Vpol=JFF[4][2];
        fc=0;
        for(k=1,#vc,
                fc=fc+vc[k]*Vpol[k];);
        return(fc);
}



EquivalentW0(JFF)=
{
        my(g,W0p);
        g=RandomPolCol(JFF,0);
        W0p=MulPointWiseFF(JFF,g,getW0(JFF));
        return(W0p);
}

FirstNonZero(liste)=
{
	my(l,flag);
	l=0;
	flag=1;
	for(k=1,#liste,
		if(flag&&(liste[k]!=0),
			l=k;
			flag=0;
			);
		);
	return(l);
}

CoefBasetoWFF(JFF,Coef,Base)=
{
	my(k0,W);
	k0=FirstNonZero(Coef);
	W=MulInt(JFF,Base[k0],Coef[k0]);
	for(k=k0+1,#Coef,
		if(Coef[k]!=0,
			W=SumJ(JFF,W,MulInt(JFF,Base[k],Coef[k]));
			);
		);
	return(W);
}

Div2Liste(liste)=
{
	my(res);
	res=vector(#liste,i,0);
	for(k=1,#liste,
		if(liste[k]%2!=0,
			print("error: WKt0even has no square root"););
		res[k]=liste[k]/2;
		);
	return(res);
}

FirstCol(liste)=
{
	return(vector(#liste,i,liste[i][1]));
}

PoltoColT(J,pol)=
{
	return(vector(poldegree(getTNicolas(J)),i,polcoef(pol,i-1))~);
}

ColtoPolT(J,c)=
{
	my(pol);
	pol=0;
	T=getTNicolas(J);
	t=variable(T);
	for(k=0,poldegree(T)-1,
		pol=pol+ c[k+1]*t^k;
		);
	return(pol);
}

\\FnDivTM(CCurve,pol)={return(FnDiv(CCurve,pol));}
FnDivTM(CCurve,pol)={return(crvfndiv(CCurve,pol));}

RiemannRochTM(CCurve,Dinit)=
{
	my(D);
	D=Dinit;
	if(type(D)=="t_VEC",D=Mat(D));
	return(crvRR(CCurve,D));
}

\\degP0Pol(pol,BranchP0)={return(BranchValuation(pol,BranchP0[1],x,y));}
degP0Pol(CCurve,pol,P)={return( crvfnord(CCurve,pol,P)  );}

degP0Col(JFF,c,CCurve,P)={return(degP0Pol(CCurve,ColtoPolFFW2(JFF,c),P));}

AlphaMap(J,W)={return(PicEval(J,W));}


FrobQq(J,elt)=
{
        my(Frobmat,taillevec,vecelt,vecFrob,eltFrob,t);
        Frobmat=J[9];
	taillevec=poldegree(getTNicolas(J));
        vecelt=vector(taillevec,i,polcoef(elt,i-1))~;\\print("vecelt = ",vecelt);
        vecFrob=Frobmat*vecelt;\\print("vecFrob = ",vecFrob);
        eltFrob=0;
        t=variable(getTNicolas(J));
        for(k=1,#vecFrob,
                eltFrob=eltFrob+(vecFrob[k]%(getpeNicolas(J)))*t^(k-1);
                );
	\\print("eltFrob = ",eltFrob);
        return(eltFrob);
}


mulpadicMatrixV2(J,s2I)=
{
	return(FqM_FqC_mul(J[11][4],s2I,getTNicolas(J),getpeNicolas(J)));
}

ColtoCoefGAdicV2(J,s)=
{
	my(Ivec,s2I,alpha);
	Ivec=J[11][5];
        s2I=vector(#Ivec,i,s[Ivec[i]])~;
        alpha=mulpadicMatrixV2(J,s2I);
	return(alpha);
}

ColtoPolGAdicV2(J,c)=
{
	my(basis,coefbasis);
	basis=J[4][2];
        coefbasis=ColtoCoefGAdicV2(J,c);
	pol=0;
	for(k=1,#coefbasis,
		pol=pol + coefbasis[k]*basis[k];
		);
	return(pol);        
}


ColtoPolCoefGAdicLoop(J,W)=
{
	my(t,kernel,flag,v,pivot);
	t=variable(W);
        for(i=1,#W[,1],for(j=1,#W,
                if(type(W[i,j])!="t_POL",W[i,j]=W[i,j] + 0*t)));
        kernel=matkerpadic(W,getTNicolas(J),getpMak(J)^geteJ(J),getpMak(J),geteJ(J));
        flag=1;
        \\print("ColtoPolCoefGAdicH : we don't assess inversibility of the pivot, only non-nullity");
	if(#kernel==0,print("ERROR: ColtoPolCoefClean: s is not in VKs"));
        for(k=1,#kernel,
                if(flag,
                        v=kernel[,k];
                        if(Mod(v[#v],getpMak(J))!=0,
                                flag=0;
                                pivot= - v[#v];
                                v=vector(#v-1,i,divpadicElt(J,v[i],pivot));
                                );
                        );
                );
        return([v,flag]);
}

ColtoPolCoefGAdicH(J,WH,c)=
\\return the polynomial of Vec(H) represented by c. We assume there exists a pol f in Vec(H) such that evaluated in Z it gives c
{
	my(W);
	W=concat(WH,c);
	return(ColtoPolCoefGAdicLoop(J,W));			
}

MatrixRemoveLign(M,i)=
{
	my(m);
	m=M;
	if(i==1,
		m=m[2..#m[,1],];);
	if(i==#M[,1],
		m=m[1..(i-1),];);
	if((i!=1)&&(i!=#M[,1]),
		m=concat(m[1..(i-1),]~,m[(i+1)..#m[,1],]~)~;);
	return(m);
}

ColtoPolCoefGAdicHWrongIndex(J,WH,c,wrongindexInit)=
\\same as ColtoPolCoefGAdich but we allow some index to be ignored
{
        my(W,wrongindex,i,res);
        W=concat(WH,c);
	wrongindex=Set(wrongindexInit);
	for(k=1,#wrongindex,
		i=wrongindex[k];
		i=i-(k-1);\\because removing ligns has changed indexes
		W=MatrixRemoveLign(W,i);
		);
	res=ColtoPolCoefGAdicLoop(J,W);
        return(res);
}

mulpadicElt(J,a,b)=
{
	return(Fq_mul(a,b,getTNicolas(J),getpMak(J)^geteJ(J)));
}

divpadicElt(J,aInit,bInit)=
{
	my(a,b);
	aInit=liftall(aInit);
	bInit=liftall(bInit);
	if(type(aInit)=="t_POL",a=aInit,
		a=aInit + 0*variable(getTNicolas(J)));
	if(type(bInit)=="t_POL",b=bInit,
		b=bInit+ 0*variable(getTNicolas(J)));
	return(ZpXQ_div(a,b,getTNicolas(J),getpMak(J)^geteJ(J),getpMak(J),geteJ(J)));
}

multPointDivSmooth(P)=
{
	my(flag,res);
	flag=1;
	res=0;
	for(k=1,#P,
		if(type(P[k])=="t_POLMOD",
			flag=0;
			res=poldegree(P[k].mod);
			);
		);
	if(flag,
		res=1;
		);
	return(res);
}

multPointDivSingular(P)=
{
	return(poldegree(P[3]));
}

multPointDiv(P)=
{
	my(res);
	if((#P==2)||(#P==3),res=multPointDivSmooth(P));
	if(#P==6,res=multPointDivSingular(P));
	return(res);
}

mulpadicMatrix(J,M1,M2)=
{
	return( FqM_mul(M1,M2,getTNicolas(J),getpeNicolas(J)) );
}

padicinvElt(J,aInit)=
{
	my(a);
	if(type(aInit)=="t_POL",a=aInit,
                a=aInit + 0*variable(getTNicolas(J)));
	return( ZpXQ_inv(a,getTNicolas(J),getpMak(J),geteJ(J)) );
}


InitFlagTF(listeThetaFraction)=
{
	my(flag);
	flag=1;
	for(k=1,#listeThetaFraction,
		if(listeThetaFraction[k][1]==0,flag=0));
	if(flag,print("GOOD Theta Fraction"),
		print("\n","ERROR: Theta Fraction","\n"));
	return(flag);
}

InitlistTF(listeThetaFraction,flagTC,flagThetaFraction)=
{
	my(listTF,listTFEven);
	if(flagTC==0,if(flagThetaFraction,listTF=vector(#listeThetaFraction,i,listeThetaFraction[i][2]),listTF=0),
		if(flagThetaFraction==0,listTF=0,
			listTF=vector(#listeThetaFraction[1],i,listeThetaFraction[1][i][2]);
			listTFEven=vector(#listeThetaFraction[2],i,listeThetaFraction[2][i][2]);
			listTF=[listTF,listTFEven];
			);
		);
	return(listTF);
}


CheckPointCurve(CCurve,Point)={return(crvcheckpt(CCurve,Point));}


CheckflagAlpha(listeRpzGx,flagTC)=
\\if flagAlphaGlobal=0, there is at least one flagAlpha which is 0
{
	my(flagAlphaGlobal,lGx);
	flagAlphaGlobal=1;
	if(flagTC,lGx=concat(listeRpzGx),
		lGx=listeRpzGx);
	for(k=1,#lGx,
		if(lGx[k][1]==0,
			flagAlphaGlobal=0));
	return(flagAlphaGlobal);
}

MCIntervertxy(PolCurveI,listPointI,P0I,P1I,P2I)=
{
	my(PolCurve,P0,P1,P2,listPoint,P);
	print("\n","Modular Curve: we intervert x and y","\n");
	PolCurve=substvec(PolCurveI,[x,y],[y,x]);
	P0=[P0I[2],P0I[1],P0I[3]];
	P1=P1I;
	if(P1I!=0,P1=[P1I[2],P1I[1],P1I[3]]);
	P2=P2I;
	if(P2I!=0,P2=[P2I[2],P2I[1],P3I[3]]);
	listPoint=vector(#listPointI,i,0);
	for(k=1,#listPointI,
		P=listPointI[k];
		P=[P[2],P[1],P[3]];
		listPoint[k]=P;
		);
	return([PolCurve,listPoint,P0,P1,P2]);
}


ThresholdOrbitCleanLoop(listsmallindex,listeFieldArchiOTCI,listeOrbiteArchiAlphaI,correspondanceI,thresholdOrbit)=
{
	my(listeFieldArchiOTC,listeOrbiteArchiAlpha,correspondance);
	print("\n","We remove the orbits of size higher than ",thresholdOrbit,"\n");
	listeFieldArchiOTC=vector(#listsmallindex,i,listeFieldArchiOTCI[listsmallindex[i]]);
        listeOrbiteArchiAlpha=vector(#listsmallindex,i,listeOrbiteArchiAlphaI[listsmallindex[i]]);
        correspondance=vector(#listsmallindex,i,correspondanceI[listsmallindex[i]]);
	\\listeRpzTFECoefBasis=vector(#listsmallindex,i,listeRpzTFECoefBasisI[listsmallindex[i]]);
	\\listeRpzGx=vector(#listsmallindex,i,listeRpzGxI[listsmallindex[i]]);
	\\listeNum=vector(#listsmallindex,i,listeNumI[listsmallindex[i]]);
	\\listeRpzLambdaGAdic=vector(#listsmallindex,i,listeRpzLambdaGAdicI[listsmallindex[i]]);
	print("listsmallindex = ",listsmallindex);
	return([listeFieldArchiOTC,listeOrbiteArchiAlpha,correspondance]);
}

ThresholdOrbitClean(listsmallindex,listeFieldArchiOTCI,listeOrbiteArchiAlphaI,correspondanceI,thresholdOrbit,flagTC)=
{
	my(listeFieldArchiOTC,listeOrbiteArchiAlpha,correspondance,listeFieldArchiOTCOdd,listeFieldArchiOTCEven,listeOrbiteArchiAlphaOdd,listeOrbiteArchiAlphaEven,correspondanceOdd,correspondanceEven);
	if(flagTC==0,[listeFieldArchiOTC,listeOrbiteArchiAlpha,correspondance]=ThresholdOrbitCleanLoop(listsmallindex,listeFieldArchiOTCI,listeOrbiteArchiAlphaI,correspondanceI,thresholdOrbit),
		[listeFieldArchiOTCOdd,listeOrbiteArchiAlphaOdd,correspondanceOdd]=ThresholdOrbitClean(listsmallindex[1],listeFieldArchiOTCI[1],listeOrbiteArchiAlphaI[1],correspondanceI[1],thresholdOrbit);
		[listeFieldArchiOTCEven,listeOrbiteArchiAlphaEven,correspondanceEven]=ThresholdOrbitClean(listsmallindex[2],listeFieldArchiOTCI[2],listeOrbiteArchiAlphaI[2],correspondanceI[2],thresholdOrbit);
		listeFieldArchiOTC=[listeFieldArchiOTCOdd,listeFieldArchiOTCEven];
		listeOrbiteArchiAlpha=[listeOrbiteArchiAlphaOdd,listeOrbiteArchiAlphaEven];
		correspondance=[correspondanceOdd,correspondanceEven];
		);
	return([listeFieldArchiOTC,listeOrbiteArchiAlpha,correspondance]);
}



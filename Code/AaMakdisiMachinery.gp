install("DivAdd","GGLGGGLL")
install("DivAdd0","GGLGGGLL")
install("DivSub","GGGLGGLGL")
install("AddFlipChain","GL")
\\install("PicIsZero","lGG") \\min
\\install("PicLiftTors","GGGLL")
\\install("PicEq","lGGG") \\min
\\install("PicFrob","GG")  \\min
install("mateqnpadic","GGGGL")
\\ \r AaTechnique.gp


\\print("SumJ: retirer le second AddFlip?? car -P=P"); On ne manipule pas toujours des points de torsion
\\Code for running Makdisi Machinery

DivAddTM(J,Wa,Wb,dim)=
{
	my(MaxAttempt,a,DA);
	MaxAttempt=50;
	a=0;
	DA=DivAdd(Wa,Wb,dim,getTNicolas(J),getpeNicolas(J),getpMak(J),a,MaxAttempt);
	if((a>=MaxAttempt)||(#DA[1,]<dim),
		\\QUEL EST LE POINTEUR NULL EN PARI 
		print("WARNING: DivAdd failed");
		a=a+1;
		DA=DivAdd0(Wa,Wb,dim,getTNicolas(J),getpeNicolas(J),getpMak(J),a,MaxAttempt););
	return(DA);
}



DivSubTM(J,Wa,Wb,dim,KV)=
\\WORKS ONLY IF Wb REPRESENTS A BASE-POINT FREE DIVISOR
\\Result has to be in Ker(KV)
{
	my(nIGS,divsub);
	nIGS=2;
	divsub=DivSub(Wb,Wa,KV,dim,getTNicolas(J),getpMak(J),geteJ(J),getpeNicolas(J),nIGS);
	if(divsub==0,
		print("ERROR DivSub");
		print("requested dim =",dim," actual dim=",#DivSubDimFF(RedJ(J,1),RedPointJ(J,Wa,1),RedPointJ(J,Wb,1)));
		);
	return(divsub);
}

DivSubTMIGS(J,Wa,Wb,dim,KV,n)={return(DivSub(Wb,Wa,KV,dim,getTNicolas(J),getpMak(J),geteJ(J),getpeNicolas(J),n));}

DivSubDimFF(JFF,Wa,Wb)=
{
	my(Ka,listeK,bj,Kj,K,S);
	Ka=WtoKFF(JFF,Wa);
	listeK=vector(#Wb[1,]+1,i,0);
	listeK[1]=Ka;
	for(j=1,#Wb[1,],
		bj=Wb[,j];
		Kj=Ka;
		Kj=matrix(#Kj[,1],#Kj[1,],k,i,MulFF(JFF,Kj[k,i],bj[i]));
		listeK[j+1]=Kj;
		);
	K=listeK[1];
	for(j=2,#listeK,
		K=concat(K~,listeK[j]~)~;);		
	S=KtoWFF(JFF,K);
	return(S);
}

DivSubDimGAdic(J,Wa,Wb)=
{
        my(Ka,listeK,bj,Kj,K,S);
        Ka=WtoKGAdic(J,Wa);
        listeK=vector(#Wb[1,]+1,i,0);
        listeK[1]=Ka;
        for(j=1,#Wb[1,],
                bj=Wb[,j];
                Kj=Ka;
                \\Kj=matrix(#Kj[,1],#Kj[1,],k,i,MulFF(JFF,Kj[k,i],bj[i]));
                Kj=matrix(#Kj[,1],#Kj[1,],k,i,mulpadicElt(J,Kj[k,i],bj[i]));
		listeK[j+1]=Kj;
                );
        K=listeK[1];
        for(j=2,#listeK,
                K=concat(K~,listeK[j]~)~;);
        S=KtoWGAdic(J,K);
        return(S);
}
	

DivSubDimKFF(JFF,Wa,Wb,KV)=
\\Rather than looking for the resulting space in Ka as usual we look for the resulting space in KV
{
        my(Ka,listeK,bj,Kj,K,S);
        Ka=WtoKFF(JFF,Wa);
        listeK=vector(#Wb[1,]+1,i,0);
        listeK[1]=KV;\\here KV
        for(j=1,#Wb[1,],
                bj=Wb[,j];
                Kj=Ka;\\BUT HERE: still Ka
                Kj=matrix(#Kj[,1],#Kj[1,],k,i,MulFF(JFF,Kj[k,i],bj[i]));
                listeK[j+1]=Kj;
                );
        K=listeK[1];
        for(j=2,#listeK,
                K=concat(K~,listeK[j]~)~;);
        S=KtoWFF(JFF,K);
        return(S);
}

DivSubTheta(J,Num,Gx,VK)=
{
	my(M,kernel,flag,v,pivot);
	M=mulPointWise(J,Gx[,1],VK);
	M=concat(M,Num);
	kernel=matkerpadic(M,getTNicolas(J),getpMak(J)^geteJ(J),getpMak(J),geteJ(J));print("DivSubTheta : dim kernel = ",#kernel);
	flag=1;
        print("DivSubTheta: matkerpadic : on ne teste pas inversibilité mais seulement la non-nullité");
        for(k=1,#kernel,
                if(flag,
                        v=kernel[,k];
                        if(v[#v]!=0,
                                flag=0;
                                pivot= - v[#v];
                                v=vector(#v-1,i,divpadicElt(J,v[i],pivot));
                                );
                        );
                );
        return([v,flag]);
}


AddFlipChainTM(n)=
{
	return(AddFlipChain(n,1));
}

AddFlipTM(J,Wa,Wb)=
\\SINON: PicChord avec flag=1
{
	my(W0,W2,K2,K3,d0,g,W12,W12p,c,W123,W3);
	W0=getW0(J);
	W2=getWn(J,2);
	K3=getKn(J,3);
	K2=getKn(J,2);
	d0=getd0(J);
	g=getg(J);
	W12=DivAddTM(J,Wa,Wb,2*d0-g+1);
	W12p=DivSubTM(J,W12,W0,d0-g+1,K3);
	c=W12p[,1];\\c=randomNNVector(W12p);
	W123=mulPointWise(J,c,W2);
	W3=DivSubTM(J,W123,W12p,d0-g+1,K2);
	return(W3);
}

AddFlipRandom(J,Wa,Wb)=
\\SINON: PicChord avec flag=1
{
        my(W0,W2,K2,K3,d0,g,W12,W12p,c,W123,W3);
        W0=getW0(J);
        W2=getWn(J,2);
        K3=getKn(J,3);
        K2=getKn(J,2);
        d0=getd0(J);
        g=getg(J);
        W12=DivAddTM(J,Wa,Wb,2*d0-g+1);
        W12p=DivSubTM(J,W12,W0,d0-g+1,K3);
	c=RandomCol(W12p);
        W123=mulPointWise(J,c,W2);
        W3=DivSubTM(J,W123,W12p,d0-g+1,K2);
        return(W3);
}

	
SumJ(J,Wa,Wb)=
{
	my(Wp,W0,W);
	Wp=AddFlipTM(J,Wa,Wb);
	W0=getW0(J);
	W=AddFlipTM(J,Wp,W0);
	return(W);
}

SumJTorsion(J,Wa,Wb)={return(AddFlipTM(J,Wa,Wb))};

SumJRandom(J,Wa,Wb)=
{
        my(Wp,W0,W);
        Wp=AddFlipRandom(J,Wa,Wb);
        W0=getW0(J);
        W=AddFlipRandom(J,Wp,W0);
        return(W);
}


MulInt(J,W,n)=
\\return [n]W
{
	my(chain,listeW,res);
	chain=AddFlipChainTM(n);
	listeW=vector(#chain+1,i,0);
	listeW[1]=getW0(J);
	listeW[2]=W;
	for(k=2,#chain,
		listeW[k+1]=AddFlipTM(J, listeW[chain[k][2][1]+1] , listeW[chain[k][2][2]+1] ););
	res=listeW[#listeW];
	return(res);
}
	

CheckZero(J,W)={return(piciszero(J,W));}
\\CheckZero(J,W)={return(PicIsZero(J,W));}
/*CheckZero(J,W)=
\\check if a point W in J rpz the origin
{
	my(W0,KW,WZ,bool);
	W0=getW0(J);
	KW=WtoK(J,W);
	WZ=DivSubDim(J,W,W0,KW);
	bool=length(WZ);
	return(bool);
}*/

CheckTorsion(J,W)={return(CheckZero(J,SumJ(J,W,W)));}	

CheckEqual(J,Wa,Wb)={return(piceq(J,Wa,Wb));}
\\CheckEqual(J,Wa,Wb)={return(PicEq(J,Wa,Wb));}
/*
CheckEgalFF(JFF,Wa,Wb)={return(#DivSubDimFF(JFF,Wa,Wb));}\\MARCHE PAS
*/

/*
WtoKGAdic(J,Wmat)=
\\EquationMatrix
\\returns the equation matrix of the space represented by W
{
\\my(listeW,listeK,K);listeW=vector(poldegree(getTNicolas(J)),i,polcoefMatrix(W,i));listeK=vector(#listeW,i,matker(W~)~);K=listeK[1];for(k=2,#listeK,K=concatEquationMatrix(K,listeK[k]);K=matimage(K~)~;);return(K);
\\rpzQqtoJMatrix(matker(rpzJtoQqMatrix(J,W~)))~
	my(W,K);
	W=rpzJtoModMatrix(J,Wmat);
	K=matker(W~)~;
	K=rpzModtoJMatrix(K);
	K=rpzQqtoJMatrix(rpzJtoQqMatrix(J,K));
	if(CheckWinK(J,Wmat,K)==0,
		print("\n","WARNING: COMPUTATION OF THE EQUATION MATRIX FAILED","\n"););
	return(K);
}*/
\\WtoKSimple(J,W)={return(rpzQqtoJMatrix(matker(rpzJtoQqMatrix(J,W~)))~);}

WtoKGAdic(J,W)={return(mateqnpadic(W,getTNicolas(J),getpeNicolas(J),getpMak(J),geteJ(J)));}

/*
KtoWGAdic(J,K)=
{
	my(W,res);
	W=matker(rpzJtoModMatrix(J,K));
	if(length(W)==0,
		res=[;];);
	if(length(W)!=0,
		W=rpzModtoJMatrix(W);
		W=rpzQqtoJMatrix(rpzJtoQqMatrix(J,W));
		if(CheckWinK(J,W,K)==0,
               		print("\n","WARNING: COMPUTATION OF W FROM THE EQUATION MATRIX FAILED","\n"););
		res=W;
		);
	return( res );
}
*/
KtoWGAdic(J,K)={return(matkerpadic(K,getTNicolas(J),getpeNicolas(J),getpMak(J),geteJ(J)));}


WtoKFF(JFF,Wmat)=
{
	my(W,K);
	W=reduceFFMatrix(JFF,Wmat);
	K=matker(W~)~;
	K=liftFFMatrix(K);
	return(K);		
}

KtoWFF(JFF,Kmat)=
{
	my(K,W);
	K=reduceFFMatrix(JFF,Kmat);
	W=matker(K);
	if(W!=[;],
		W=liftFFMatrix(W);
		);
	return(W);
}

\\LiftTorsion(J,W)={return(PicLift(J,W,1,2,0));}
LiftTorsion(J,W)={return(piclift(J,W,1,2,0));}

FrobMakdisi(J,W)={return(picfrob(J,W));}	





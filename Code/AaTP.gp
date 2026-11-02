\\ \r AaMakdisiMachinery.gp
\\ \r AaTechnique.gp

\\We interpolate the theta polynomials




InitArchidataCore(J,ChiOTC,listeOrbiteAlpha,flagTC,thresholdOrbit)=
\\Identification of the defnition field of each OTC
{
        my(chi,listeField,listeFieldArchiOTC,zv,DefFieldVec,listeOrbiteArchiAlpha,correspondance,alpha,OrbiteGAdicalpha,flag,d,deg,flagTest,alphatest,a,b,listsmallindex,flagThreshold);
        \\if(flagTC,print("	EVEN TC computation"));
	flagThreshold=0;\\if flag=1, we have remove some orbit(s)
	chi=ChiOTC;
        listeField=factor(chi)[,1];
        listeFieldArchiOTC=vector(#listeField,i,0);\\definition Fields of the Alpha
	\\zv='z;
	zv=varlower("w",y);
	listsmallindex=[];\\gather the index of the orbits whose size is smaller than the threshold. There are the ones we keep.
	for(i=1,#listeFieldArchiOTC,
		if(poldegree(listeField[i])>thresholdOrbit,flagThreshold=1,
			listsmallindex=concat(listsmallindex,[i]);
			DefFieldVec=nfinit(subst(listeField[i],variable(listeField[i]),zv),3);
			if(#DefFieldVec!=2,\\the polynomial given was monomial, so nf is a "normal" nf
				DefFieldVec=[DefFieldVec,Mod(zv,DefFieldVec.pol)];
				);
			listeFieldArchiOTC[i]=DefFieldVec;
			);
		);
	listeOrbiteArchiAlpha=vector(#listeField,i,[]);
	correspondance=vector(#listeField,i,[]);
	\\listeOrbiteArchiAlpha[k] = (list of list) one orbite over Q which is decomposed in list of elements in the same Frob/GAdic orbit
	\\listeOrbiteArchiAlpha[k][i] = (list) subset of conjuagtes in the same Frob orbit
	\\listeOrbiteArchiAlpha[k][i][j] =  element rpz alpha as a polynomial in Qp/p^e[t] 
        for(k=1,#listeOrbiteAlpha,
                alpha=Mod(listeOrbiteAlpha[k][1]*(1+O( getpMak(J)^geteJ(J) )), getTNicolas(J) );
		OrbiteGAdicalpha=listeOrbiteAlpha[k];
                flag=1;
                for(i=1,#listeField,
                        if(flag,
                                alphatest=alpha*(1+O(getpMak(J)^(emax-20)));
				if(liftall(substpol(listeField[i],x,alphatest))==0,
                                        flag=0;
					correspondance[i]=concat(correspondance[i],[k]);
                                        listeOrbiteArchiAlpha[i]=concat(listeOrbiteArchiAlpha[i],[OrbiteGAdicalpha]);
					);
                                );
                        );
                if(flag,
                        print("WARNING: InitlisteDefFieldRpzOrbiteOTC, identification of one of the alpha failed");
			print("k = ",k);
			print("indexes that we aim to keep: ",listsmallindex);
			print("alpha = ",alpha*(1+O(getpMak(J)^2)));
			for(i=1,#listField,
				b=alpha*(1+O(getpMak(J)^2));
				a=lift(substpol(listeField[i],x,b));
				print(a);
				print("IsZero? ",a==0)	
				);
			);
                );
	flagTest=1;
	for(k=1,#listeField,
		deg=poldegree(listeField[k]);
		d=0;
		for(i=1,#listeOrbiteArchiAlpha[k],
			d=d+#listeOrbiteArchiAlpha[k][i];
			);
		if(d!=deg,
			flagTest=0;
			);
		);
	if(default(debug)>1,print("Test listeOrbiteArchiAlpha : ",flagTest));

        return([listeFieldArchiOTC,listeOrbiteArchiAlpha,correspondance,listsmallindex,flagThreshold]);
}

InitArchidata(J,ChiOTC,listeOrbiteAlpha,flagTC,thresholdOrbit)=
{
	my(listeFieldArchiOTC,listeOrbiteArchiAlpha,correspondance,listeFieldArchiOTCEven,listeOrbiteArchiAlphaEven,correspondanceEven,listsmallindex,listsmallindexEven,flagT,fT);
	if(default(debug)>1,print("\n","beginning InitArchidata"));
	gettime();
	if(flagTC,[listeFieldArchiOTC,listeOrbiteArchiAlpha,correspondance,listsmallindex,flagT] = InitArchidataCore(J,ChiOTC[1],listeOrbiteAlpha[1],0,thresholdOrbit),
		[listeFieldArchiOTC,listeOrbiteArchiAlpha,correspondance,listsmallindex,flagT] = InitArchidataCore(J,ChiOTC,listeOrbiteAlpha,0,thresholdOrbit));

	if(flagTC,
		[listeFieldArchiOTCEven,listeOrbiteArchiAlphaEven,correspondanceEven,listsmallindexEven,fT] = InitArchidataCore(J,ChiOTC[2],listeOrbiteAlpha[2],1,thresholdOrbit);
		listeFieldArchiOTC=[listeFieldArchiOTC,listeFieldArchiOTCEven];
		listeOrbiteArchiAlpha=[listeOrbiteArchiAlpha,listeOrbiteArchiAlphaEven];
		correspondance=[correspondance,correspondanceEven];
		listsmallindex=[listsmallindex,listsmallindexEven];
		if(fT==1,flagT=1);
		);
	if(default(debug)>1,print("end InitArchidata. time computation: ",strtime(gettime())));
	return([listeFieldArchiOTC,listeOrbiteArchiAlpha,correspondance,listsmallindex,flagT]);
}

InitW2plusP0(J,CCurve,P0)=
\\H°(2D0+P0)
{
	my(D,H2plusP0,W2plusP0);
	print("if alpha is well defined everywhere, please skip this step by setting flagAlphaGlobal=1");
	D=matrix(1,2,i,j,0);
	D[1,1]=P0;
	D[1,2]=2*getd0(J)+1;
	H2plusP0=RiemannRochTM(CCurve,D);
	W2plusP0=PoltoMakdisi(J,H2plusP0);
	return(W2plusP0);
}

\\InitlisteThetaPolynomial
export(DivSubDimFF)

InitlisteRpzGxCore(J,listeRpzOrbiteOTC,parity,W2plusP0)=
\\ (Gx) = D + E1 + Ex - 2D0
\\ Gx is (generically) the unique function of H°(2D0 - D - E1)
{
	my(WE1,W2,KV,d0,listeRpzGx,WD,Wp,W,Gx,flagAlpha);
	if(parity==0,print("	EVEN TC computation"),print("	ODD TC computation"));
	WE1=getWE(J)[1][1];\\J[14][1][1];\\H°(2D0 - E1)
	W2=getWn(J,2);
	KV=getKn(J,2);
	d0=getd0(J);
	listeRpzGx=vector(#listeRpzOrbiteOTC,i,0);

/*
	for(k=1,#listeRpzGx,
		print("InitlisteRpzGX: new representant (",k,"/",#listeRpzGx,")"); 
		WD=listeRpzOrbiteOTC[k];
		Wp=DivAddTM(J,WD,WE1,2*d0+1);\\H°(4D0 - D - E1)
		W=DivSubTM(J,Wp,W2,1,KV);\\H°(2D0 - D - E1)       = H°(2D - D - E1 - Ex)
		if(W==0,print("\n"," error: try with another WE1/WE2","\n"));
		Gx=W[,1];
		listeRpzGx[k]=Gx;
		);
	\\A=listeRpzGx;
*/
	if(default(debug)>1,print("InitlisteRpzGx"));

	parfor(k=1,#listeRpzGx,
                print("new representative (",k,"/",#listeRpzGx,")");
		flagAlpha=1;\\if flagAlpha=0, then we had to impose more constraints on P0 to define the Alpha map
                WD=listeRpzOrbiteOTC[k];
                Wp=DivAddTM(J,WD,WE1,2*d0+1);\\H°(4D0 - D - E1)
                W=DivSubTM(J,Wp,W2,1,KV);\\H°(2D0 - D - E1)       = H°(2D - D - E1 - Ex)
                if(W==0,

			print("\n"," alpha not defined on this point: we impose one more condition on P0","\n");
			flagAlpha=0;
			W=DivSubTM(J,Wp,W2plusP0,1,KV);\\H°(2D0 -P0 - D - E1)
			if(W==0,print("\n","we still have an error","\n"));
			if(W!=0,print("dim Gx == 1 ?",#W));
			);

                Gx=W[,1];
		[flagAlpha,Gx],
			r,listeRpzGx[k]=r);

	\\B=listeRpzGx;print("A==B ??",A==B);

	print("end");
	return(listeRpzGx);
}

InitlisteRpzGx(J,listeRpzOrbiteOTC,flagTC,W2plusP0)=
\\ (Gx) = D + E1 + Ex - 2D0
\\ Gx is (generically) the unique function of H°(2D0 - D - E1)
{
	my(listeRpzGx,listeRpzGxEven);
	if(flagTC,listeRpzGx=InitlisteRpzGxCore(J,listeRpzOrbiteOTC[1],1,W2plusP0),
		listeRpzGx=InitlisteRpzGxCore(J,listeRpzOrbiteOTC,1,W2plusP0));

	if(flagTC,
		listeRpzGxEven=InitlisteRpzGxCore(J,listeRpzOrbiteOTC[2],0,W2plusP0);
		listeRpzGx=[listeRpzGx,listeRpzGxEven];
		);
	return(listeRpzGx);
}

InitHK0tRR(J,Omega,CCurve,P0)=
\\return a basis of polynomial of H°(K0t) with K0t = K0 + (d0-(2g-2))P0
\\use RiemannRoch computation
{
	my(d,g,counter,c,f,Df,index,DK0,flag,P0p,D,HK0t,VK0t,dif,flag2);
        d=0;
        g=getg(J);
	counter=0;
        while(d!=(2*g-2),
		counter=counter+1;
		if(default(debug)>2,print("HK0t RR: counter = ",counter," d = ",d));
                d=0;
                c=vector(#Omega,i,1+random(30))~;
                f=Omega*c;
                Df=FnDivTM(CCurve,f);\\print("Df = ",Df);
                index=1;
                while(Df[index,2]>=0,
                        index=index+1;
                        );
                DK0=Df[index..#Df[,1],];\\print("DK0 = ",DK0);
		for(k=1,#DK0[,2],
			DK0[k,2]=-DK0[k,2];
			\\print("DK0[k,2]",DK0[k,2]);
			);
		d=crvdivdeg(CCurve,DK0);
                \\for(k=1,#DK0[,1],
                \\        d=d+DK0[k,2]*multPointDiv(DK0[k,1]);
                \\        );
                );
        \\print("d=",d);
        flag=0;
        P0p=P0;
        if(#P0==2,
                P0p=[P0[1],P0[2],1];);
	\\print("DK0 =",DK0);
	
	dif=getd0(J) - (2*getg(J) -2);
	flag2=1;
	D=DK0;\\print("P0p = ",P0p);print("DK0 = ",DK0[1,1]);
	for(k=1,#DK0[,1],
		if(CheckEgalitePointCoord(P0p,DK0[k,1]),
			flag2=0;
			D[k,2]=DK0[k,2]+dif;
			);
		);
	if(flag2,
        	D=concat([P0p,dif]~,DK0~)~;
		);
	if(default(debug)>2,print("K0t =",D));
	
        HK0t=RiemannRochTM(CCurve,D);\\print("HK0t = ",HK0t);
        VK0t=PoltoMakdisi(J,HK0t);\\evaluation matrix of the space H°(K0t). Do not confuse with WK0t
        return([HK0t,VK0t,D]);
}

/*
InitHK0tD0(J,HK0t)=
\\return a basis of polynomials for H°(D0 + K0t)
{
	my();
	HD0=J[4][1];
	\\nIGS=2;
	dim=2*getd0(J) -getg(J) +1;
	d=0;
	while(d<dim,
		f=RandomPolH(HK0t);
		g=RandomPolH(HK0t);
		Hf=vector(#HD0,i,f*HD0[i]);
		Hg=vector(#HD0,i,g*HD0[i]);
		W=concat(PoltoMakdisi(J,Hf),PoltoMakdisi(J,Hg));
		W=\\WE NEED A MATIMAGE IN J(Zq/p^e)
		);
}
*/

InitHK0tD0RR(J,DK0t,CCurve,P0)=
\\return a basis of polynomials for H°(D0 + K0t)
\\use RiemannRoch
{
	my(P0p,flag2,D,HK0tD0);
	P0p=P0;
        if(#P0==2,
                P0p=[P0[1],P0[2],1];);
	flag2=1;
        D=DK0t;
        for(k=1,#DK0t[,1],
                if(CheckEgalitePointCoord(P0p,DK0t[k,1]),
                        flag2=0;
                        D[k,2]=DK0t[k,2]+getd0(J);
                        );
                );
        if(flag2,
                D=concat([P0p,getd0(J)]~,DK0t~)~;
                );
	if(default(debug)>2,print("DK0tD0 = ",D));
	
	HK0tD0=RiemannRochTM(CCurve,D);
	return(HK0tD0);
}

InitK_K0tD0(J,HK0tD0)=
\\K the equation matrix for VK0tD0 with represents H°(D0 + K0t)
{
	my(VK0tD0,K_K0tD0);
	VK0tD0=PoltoMakdisi(J,HK0tD0);
	K_K0tD0=WtoKGAdic(J,VK0tD0);
	return([K_K0tD0,VK0tD0]);
}








\\ Num
ThetaPolynomialNum(J,WOTC,VK0t,K_K0tD0)=
\\ (Num) = 2D - D0 - K0t       ou [D-D0]=OTC tel qu'il est donné
{
	my(d0,g,W0,Wp,Wpp,Wppp,W,Num);
	\\print("InitlisteNum: new loop");
	d0=getd0(J);
	g=getg(J);
	W0=getW0(J);
	Wp=DivAddTM(J,WOTC,WOTC,2*d0 -g +1);\\H°(4D0-2D)
	Wpp=DivSubTM(J,Wp,W0,d0 -g +1,getKn(J,3));\\H°(3D0-2D)
	\\print("Wpp done");
	Wppp=DivAddTM(J,Wpp,VK0t,2*d0 -g +1);\\H°(3D0+K0t-2D)? the number of evaluation points should be enough to represent a unique point in H°(3D0+K0t). Maybe this matrice represents another point in H°(5D0) but we don't care since we don't need to distinguish.
	\\print("Wppp done");
	if(Wppp==0,print("\n"," error: try with another WE1/WE2","\n"));
	W=DivSubTM(J,Wppp,getWn(J,2),1,K_K0tD0);\\H°(D0+K0t-2D)
	\\print("W done");
	if(W!=0,Num=W[,1],
		print("\n error DivSub ThetaPolynomialNum (listeNum) \n");
		Num=0;
		);
	return(Num);
}

export(ThetaPolynomialNum)

InitlisteNumCore(J,listeRpzOrbiteOTC,VK0t,K_K0tD0,parity)=
{
	my(listeNum);
	if(default(debug)>1,print("InitlisteNum"));
	if(parity==0,print("	EVEN TC computation"),print("	ODD TC computation"));

	listeNum=vector(#listeRpzOrbiteOTC,i,0);

/*
	for(i=1,#listeRpzOrbiteOTC,
		print("InitlisteNum: new loop (",i,"/",#listeNum,")");
		listeNum[i]=ThetaPolynomialNum(J,listeRpzOrbiteOTC[i],VK0t,K_K0tD0);
		);
	A=listeNum;
*/

	parfor(i=1,#listeRpzOrbiteOTC,
                print("new representative (",i,"/",#listeNum,")");
		ThetaPolynomialNum(J,listeRpzOrbiteOTC[i],VK0t,K_K0tD0),
			r,listeNum[i]=r);

	\\B=listeNum;print("A==B ??",A==B);

	print("end");
	return(listeNum);
	\\return(vector(#listeRpzOrbiteOTC,i,ThetaPolynomialNum(J,listeRpzOrbiteOTC[i],VK0t,K_K0tD0)));
}

InitdataRRspace(J,Omega,CCurve,P0)=
{
	my(HK0t,VK0t,DK0t,HK0tD0,K_K0tD0,VK0tD0,H3D0,W3,RR3D0);
	if(default(debug)>1,print("beginning RR computation"));
	gettime();
	[HK0t,VK0t,DK0t]=InitHK0tRR(J,Omega,CCurve,P0);if(default(debug)>1,print("HK0t : okay"));
        HK0tD0=InitHK0tD0RR(J,DK0t,CCurve,P0);if(default(debug)>1,print("HK0tD0 : okay"));
        [K_K0tD0,VK0tD0]=InitK_K0tD0(J,HK0tD0);if(default(debug)>1,print("K_K0tD0 : okay"));
        \\VK0tD0=PoltoMakdisi(J,HK0tD0);
        \\VK0tD0="we could compute VK0tD0 but unless this message is displaying, we don't need it";
	H3D0=RiemannRochTM(CCurve,Mat([P0,3*getd0(J)])); \\TODO: computation could be upgraded using H2
        if(default(debug)>1,print("end"));
	W3=PoltoMakdisi(J,H3D0);
        RR3D0=[H3D0,W3];
        dataRRspace=[HK0tD0,K_K0tD0,VK0tD0,HK0t,VK0t,DK0t,RR3D0];	

	print("end dataRRspace. time computation: ",strtime(gettime()));
	return(dataRRspace);
}

InitlisteNum(J,listeRpzOrbiteOTC,dataRRspace,flagTC)=
{
	my(listeNum,HK0t,VK0t,DK0t,HK0tD0,K_K0tD0,VK0tD0,listeNumEven,H3D0,W3,RR3D0);

	[HK0tD0,K_K0tD0,VK0tD0,HK0t,VK0t,DK0t,RR3D0]=dataRRspace;
	[H3D0,W3]=RR3D0;

	if(flagTC,listeNum=InitlisteNumCore(J,listeRpzOrbiteOTC[1],VK0t,K_K0tD0,1),
		listeNum=InitlisteNumCore(J,listeRpzOrbiteOTC,VK0t,K_K0tD0,1));

	if(flagTC,
		listeNumEven=InitlisteNumCore(J,listeRpzOrbiteOTC[2],VK0t,K_K0tD0,0);
		listeNum=[listeNum,listeNumEven];
		);

	return(listeNum);
}







\\Theta Polynomial
/*
InitsTheta(J,Gx,Num)=
\\since Num/Gx^2 = fBeta a ThetaPolynomial as a polynomial and that the divisor of fBeta is (fBeta)=2Beta-K0, we have that Num/Gx^2 is in H°(D0+K0t). So even if it is not generically possible to divide to element in a Makdisi framework, we should be able to do it in that situation
{
	my(sTheta,flag,gx);
	sTheta=vector(#Gx,i,0)~;
	flag=1;
	for(i=1,#Gx,
		gx=Gx[i];
		if(gx == 0, \\then Z[i] is either also a 0 of Num or a point of K0. But this second situation is supposed avoided due to the constuction of Z
			flag=0;
			print("Num et Gx ont un 0 commun dans les points d'eval: A CODER ou renvoyer le theta polynomial sous forme de fraction rationnelle");
			);
		if(flag,
			sTheta[i]=divpadicElt(J,Num[i],mulpadicElt(J,gx,gx));
			);
		);
	return([sTheta,flag]);
}
*/

InitsThetaDivSub(J,Gx,Num,VK0tD0)=
{
	my(Gxm,Numm,sTheta,flag);
	flag=0;
	Gxm=matrix(#Gx,1,i,j,mulpadicElt(J,Gx[i],Gx[i]));\\Gx^2
	Numm=matrix(#Num,1,i,j,Num[i]);
	\\sTheta=DivSubTM(J,Numm,Gxm,1,WtoKGAdic(J,VK0tD0))[,1];
	\\print("if flagAlphaGlobal = 0, it may bug here ??   BIS  ");
	[sTheta,flag]=DivSubTheta(J,Numm,Gxm,VK0tD0);
	if(flag,
		print("DivSubTheta : error: sTheta = ",sTheta);
		);
	return(sTheta);
}

InitsThetaDivSubInv(J,Gx,Num)=
{
        my(Gxm,Numm,sThetaInv);
        Gxm=matrix(#Gx,1,i,j,mulpadicElt(J,Gx[i],Gx[i]));\\Gx^2
        Numm=matrix(#Num,2,i,j,Num[i]);
	\\print("if flagAlphaGlobal = 0, it may bug here ??");
        sThetaInv=DivSubTM(J,Gxm,Numm,1,getKn(J,3))[,1];\\Gx^2/Num
        return(sThetaInv);
}

/*
CoefArchialgdep(coefGAdic,DefField)=
{
	my(t,T,coefGAdicElt,polmin,nfpol,basis,coefArchi);
	t=varlower("t",x);
	T=getTNicolas(J);
	coefGAdicElt=Mod(subst(coefGAdic,variable(coefGAdic),t),subst(T,variable(T),t));\\print("coefGAdicElt = ",coefGAdicElt);
	polmin=bestappr(algdep(coefGAdicElt,poldegree(DefField.pol)));
	nfpol=subst(DefField.pol,x,t);
        basis=vector(#poldegree(nfpol),i,Mod(subst(DefField.zk[i],x,t),nfpol));
	coefArchi=nflindep(concat([coefGAdicElt],basis));
	if(coefArchi[1]==0,
        	print("Interpolation Theta: a coefficient has not been identificated");
                coef=-1;
                );
	if(coefArchi[1]!=0,
                coef=vector(#coefArchi-1,i,-coefArchi[i+1]/coefArchi[1])~;
                );
	coefArchi=nfbasistoalg(DefField,coefArchi);
	return(coefArchi);
}
*/

/*
InterpolationCoeflindep(J,cGAdic,DefField,alpha)=
{
	my(polnf,n,iprim,T,t,basisArchi,cArchiCoef,pivot,cArchiCoefpropre,a,cArchi,cible,basis,l);
	polnf=DefField.pol;
	n=poldegree(polnf);
	iprim=alpha*(1+O(getpMak(J)^geteJ(J)));\\polynome in Zp/p^e which fairly represents the archimedian element Mod(x,DefField.pol) in Zq
	T=getTNicolas(J);	
	t=variable(T);
	iprim=subst(iprim,variable(iprim),t);
	print("InterpolationCoef : injection : ",liftall(subst(polnf,x,Mod(iprim,T))) == 0);
	if(liftall(subst(polnf,x,Mod(iprim,T)))!=0,
                print("InterpolationCoef : error injection");
                );
	basisArchi=vector(n,k,lift(Mod(iprim,T)^(k-1)));
	\\cArchiCoef=lindep(concat(basisArchi,[cGAdic*(1+O(getpMak(J)^geteJ(J)))]));print("InterpolationCoef : cArchiCoef = ",cArchiCoef);
	cArchiCoef=lindep(concat([cGAdic*(1+O(getpMak(J)^geteJ(J)))],basisArchi));print("InterpolationCoef : cArchiCoef = ",cArchiCoef);

		
	\\if we give polynomials to lindep, the returned coef are in the ground field and not integers anymore. So we do degree by degree
	\\cGadic=cGAdic*(1+O(getpMak(J)^geteJ(J)));
	\\cible=polcoef(cGAdic,0);
	\\pivot=0;
	\\basis=vector(#basisArchi,i,polcoef(basisArchi[i],0));
	\\cArchiCoef=vector(#basisArchi,i,0);\\integers
	\\for(k=0,poldegree(T)-1,
	\\	l=lindep(concat([cible],basis));
	\\	pivot=l[1];
		
	
	print("bestappr : ",bestappr(cArchiCoef));	

	\\cArchiCoef=bestappr(cArchiCoef);print("InterpolationCoef : cArchiCoef post bestappr = ",cArchiCoef);
	pivot=cArchiCoef[1];
	if(lift(pivot)==0,
		print("InterpolationCoef: error pivot = 0");
		);
	\\pivot=pivot*(1+O(getpMak(J)^geteJ(J)));
	\\cArchiCoef=vector(#cArchiCoef,i,cArchiCoef[i]/(-pivot));\\last entry is -1
	cArchiCoefpropre=vector(#cArchiCoef,i,0);
	for(i=1,#cArchiCoef,
		\\a=cArchiCoef[i]*(1+O(getpMak(J)^geteJ(J)));
		a=cArchiCoef[i];
		a=a/(-pivot);
		\\a=divpadicElt(J,lift(a)+0*t,lift(-pivot)+0*t)*(1+O(getpMak(J)^geteJ(J)));
		a=bestappr(a);
		if(a==[],
			print("WARNING: Interpolation Coef: bestappr failed or a=0");
			);
		cArchiCoefpropre[i]=a;
		);
	cArchiCoef=cArchiCoefpropre;
	print("InterpolationCoef : cArchiCoefpropre = ",cArchiCoef);
	cArchi=0;
	for(k=0,#cArchiCoef-2,
		cArchi=cArchi + cArchiCoef[k+2]*x^k;
		);
	cArchi=Mod(cArchi,polnf);print("InterpolationCoef: cArchi = ",cArchi);
	return(cArchi);
}
*/

InterpolationCoefLagrange(J,lOrbc,DefFieldVec,OrbiteArchiAlpha)=
\\ATTENTION: NE PEUT PAS MARCHER SI DefField n'est pas galoisien  OU ALORS prende toute l'orbite archimédienne
{
	my(DefField,cDefField,T,t,dT,n,X,p,e,OrbitecGAdic,Y,taille,index,prepol,poly,cp,cArchi,flag);
	[DefField,cDefField]=DefFieldVec;
	\\z=varlower("z",y);
	\\nfz=subst(DefField.pol,x,z);
	T=getTNicolas(J);
	t=variable(T);
	dT=poldegree(T);
	n=#OrbiteArchiAlpha;
	X=concat(OrbiteArchiAlpha);
	p=getpMak(J);
	e=geteJ(J);
	X=vector(#X,i,Mod(X[i]*(1+O(p^e)),T));
	Y=vector(n,i,0);
	for(i=1,n,
		taille=#OrbiteArchiAlpha[i];
		OrbitecGAdic=vector(taille,j,0);
		OrbitecGAdic[1]=lOrbc[i];
		for(j=2,taille,
                	OrbitecGAdic[j]=FrobQq(J,OrbitecGAdic[j-1]);
                	);
		OrbitecGAdic=vector(#OrbitecGAdic,j,Mod(OrbitecGAdic[j]*(1+O(p^e)),T));
		Y[i]=OrbitecGAdic;
		);
	Y=concat(Y);
	prepol=polinterpolate(X,Y,'x);\\print("InterpolationArchi: prepol = ",prepol);
	poly=0;
	flag=1;
	for(k=0,poldegree(prepol,'x),
		c=polcoef(prepol,k,'x);
		c=lift(c);
		if((c!=0)&&(vector(poldegree(c,t),j,polcoef(c,j,t))!=0),
			print("error InterpolationCoefLagrange: coef not in Qp : ",vector(poldegree(c,t),j,polcoef(c,j,t))," degree = ",poldegree(c,t));
			flag=0;
			cp=0;
			);
		c=polcoef(c,0,t);
		cp=bestappr(c);\\print(" InterpolationCoefLagrange: cp = ",cp," lift(c) = ",lift(c));
		if((c!=0)&&(cp==[]),
			cp=0;
			flag=0;
			print(" InterpolationCoefLagrange, ERROR bestappr: k= ",k);  \\lift(c));
			);
		poly=poly+cp*y^k;
		);
	prepol=bestappr(poly);\\print("InterpolationArchi: prepol post bestappr = ",prepol);
	\\cArchi=Mod(prepol,nfz);we have to take account of the case where nf is given by nfinit(non-monomial polynomial)
	cArchi=subst(prepol,y,cDefField);\\print("cArchi = ",cArchi);
	return([flag,cArchi]);
}

InterpolationTheta(J,lOrbsCoefBasis,DefFieldVec,HKs,OrbiteArchiAlpha)=
{
	my(sArchi,lOrbc,cArchi,flag,flagcoef);
	\\t=variable(getTNicolas(J));
	sArchi=0;
	flag=1;
	for(k=1,#HKs,
		\\coefGAdic=sGAdic[k]*(1+O(getpMak(J)^geteJ(J)));
		\\coefArchi=CoefArchialgdep(coefGAdic,DefField);
		lOrbc=vector(#lOrbsCoefBasis,j,lOrbsCoefBasis[j][k]);
		[flagcoef,cArchi]=InterpolationCoefLagrange(J,lOrbc,DefFieldVec,OrbiteArchiAlpha);
		if(flagcoef==0,flag=0);
		sArchi=sArchi + cArchi*HKs[k];
		);
	return([flag,sArchi]);
}
	
/*	
InterpolationColGAdic(J,s,HKs,VKs,DefField,OrbiteArchiAlpha)=
{
	my(sGAdic,flag,lpivot,pivot,sArchi);
	[sGAdic,flag]=ColtoPolCoefGAdicH(J,HKs,VKs,s);print("CoefGAdic done");
        if(flag==0,
        	\\on rigidifie
                lpivot=CoefNonZero(sGAdic);
                pivot=sGAdic[lpivot[#lpivot]];
                sGAdic=vector(#sGAdic,i,divpadicElt(J,sGAdic[i],pivot));
                sArchi=InterpolationTheta(J,sGAdic,DefField,HKs,OrbiteArchiAlpha);print("sArchi done");
                );
	if(flag,
        	print("error in InterpolationTheta");
                );
	return(sArchi);
}
*/

ColtoPolCoefClean(J,s,VKs,wrongindex)=
{
        my(sGAdic,flag,lpivot,pivot);
	if(wrongindex==0,
		[sGAdic,flag]=ColtoPolCoefGAdicH(J,VKs,s);if(default(debug)>2,print("CoefGAdic done"));
		);
	if(wrongindex!=0,
		if(default(debug)>1,print("sTheta[k] not always invertible: 'wrongindex' techniques"));
		[sGAdic,flag]=ColtoPolCoefGAdicHWrongIndex(J,VKs,s,wrongindex);if(default(debug)>2,print("CoefGAdic done"));
                );
	if(flag,flag=0,flag=1);\\everything's good if flag=1
        if(flag,  \\VKs*sGAdic = s
                \\on rigidifie
                \\lpivot=CoefNonZero(sGAdic);
		lpivot=CoefInversible(sGAdic,getpMak(J));
		if(default(debug)>2,print("pivot's position (expected unchanged in the whole Global orbit): ",lpivot[#lpivot]));
                pivot=sGAdic[lpivot[#lpivot]];\\normalement, puisque le pivot est le dernier coef non nul, sa position devrait être la même dans une même orbite archi
                sGAdic=vector(#sGAdic,i,divpadicElt(J,sGAdic[i],pivot));
                );
        if(flag==0,
                print("\n","error in InterpolationTheta (CoeftoPolCoefClean)","\n");
                );
        return([flag,sGAdic]);
}


/*
ThetaPolynomial(J,WOTC,Gx,HK0tD0,VK0t,HK0t,K_K0tD0,DefField,Num,OrbiteArchiAlpha)=
{
	my(sTheta,flags,Theta,sThetaInv);
	print("ThetaPolynomial: nouvelle boucle");
	[sTheta,flags]=InitsTheta(J,Gx,Num);print("InitsTheta done");\\Column rpz Theta in Makdisi framework
	\\if(flags==0,
		
		\\print("on recours à InitsThetaDivSub");
		\\sTheta=InitsThetaDivSub(J,Gx,Num,K_K0tD0);
		\\print("	    done");
		
		\\fraction rationnelle ou calcul du degré de Num et Gx en les 0 communs dans Z ?
		\\);
	if(flags,
		Theta=InterpolationColGAdic(J,sTheta,HK0t,VK0t,DefField,OrbiteArchiAlpha);\\pol
		);
	if(flags==0,
		\\sTheta=InitsThetaDivSub(J,Gx,Num,K_K0tD0);
		\\Theta=InterpolationColGAdic(J,sTheta,HK0t,VK0t,DefField,alpha);
		print("Theta sous forme de fraction rationnelle");
		sThetaInv=InitsThetaDivSubInv(J,Gx,Num);\\return la colonne qui rpz 1/Theta in H°(2D0)
		ThetaInv=InterpolationColGAdic(J,sThetaInv,J[4][2],getWn(J,2),DefField,OrbiteArchiAlpha);\\return le pol Archi 1/Theta
		Theta=1/ThetaInv;
		);	
	return(Theta);		
}
*/

listeThetaPolynomialCore(J,listeRpzGx,listeNum,listeOrbiteArchiAlpha,listeFieldArchiOTC,correspondance,listeTFE,failorbit,RR3D0)=
{
	my(H3D0,W3,flaginv,listeThetaPolynomial,lindex,OrbiteArchiAlpha,Field,lOrbsThetaInv,lOrbsThetaInvCoefBasis,flagInterpolation,ThetaInv,lOrbsTheta,lOrbsThetaCoefBasis,Theta);
	flaginv=1;
	listeThetaPolynomial=listeTFE;
	\\[HK0tD0,K_K0tD0,VK0tD0,HK0t,VK0t,DK0t]=dataRRspace;
	\\H2D0=J[4][2];
	\\H2D0=getW2Pol(J);
	\\W2=getWn(J,2);
	[H3D0,W3]=RR3D0;\\Becareful!! 1/f is in H°(3D0) in general and not in H°(2D0)
	for(i=1,#listeThetaPolynomial,
		if(IsInList(failorbit,i),
			print("Theta Polynomial: new Global orbit (",i,"/",#listeThetaPolynomial,")");
			lindex=correspondance[i];
			\\lOrbGx=vector(#lindex,k,listeRpzGx[lindex[k]]);
			\\lOrbNum=vector(#lindex,k,listeNum[lindex[k]]);
			OrbiteArchiAlpha=listeOrbiteArchiAlpha[i];
			Field=listeFieldArchiOTC[i];
			if(flaginv,
				print("Theta as rational fraction");
				lOrbsThetaInv=vector(#lindex,k,InitsThetaDivSubInv(J,listeRpzGx[lindex[k]],listeNum[lindex[k]]));
				lOrbsThetaInvCoefBasis=vector(#lindex,k,ColtoPolCoefClean(J,lOrbsThetaInv[k],W3,0)[2]);
				[flagInterpolation,ThetaInv]=InterpolationTheta(J,lOrbsThetaInvCoefBasis,Field,H3D0,OrbiteArchiAlpha);
				if(flagInterpolation,
					Theta=1/ThetaInv;
					);
				if(flagInterpolation==0,
					print("One Theta Fraction failed"));
		\\if(flaginv==0,
		\\	lOrbsTheta=vector(#lindex,k,InitsThetaDivSub(J,listeRpzGx[lindex[k]],listeNum[lindex[k]],VK0tD0));
		\\	lOrbsThetaCoefBasis=vector(#lindex,k,ColtoPolCoefClean(J,lOrbsTheta[k],W2));
		\\	Theta=InterpolationTheta(J,lOrbsThetaCoefBasis,Field,H2D0,OrbiteArchiAlpha);
		\\	);
				listeThetaPolynomial[i]=[flagInterpolation,Theta];
				);
			);
		);
	return(listeThetaPolynomial);
}
/*
InitlisteThetaPolynomial(J,listeRpzOrbiteOTC,listeFieldArchiOTC,listeOrbiteArchiAlpha,Omega,CCurve,P0,correspondance)=
\\ (Theta) = 2Btt - D0 - K0t ou Btt est un diviseur archimédien tq Btt-(d0 -(g-1))P0 ~ Beta un OTC
{
	my(H2D0,listeRpzGx,HK0t,VK0t,DK0t,HK0tD0,K_K0tD0,VK0tD0,listeNum,listeTFE,failorbit,listeThetaPolynomial);
	\\H2D0=J[4][2];
	H2D0=getW2Pol(J);
	listeRpzGx=InitlisteRpzGx(J,listeRpzOrbiteOTC);\\still under Makdisi representation
	print("listeRpzGx: done");
	[HK0t,VK0t,DK0t]=InitHK0tRR(J,Omega,CCurve,P0);print("HK0t : done");
	print("long ? (InitHK0tD0RR)");
	HK0tD0=InitHK0tD0RR(J,DK0t,CCurve,P0);print("HK0tD0 : done");
	K_K0tD0=InitK_K0tD0(J,HK0tD0);print("K_K0tD0 : done");
	VK0tD0=PoltoMakdisi(J,HK0tD0);
	listeNum=InitlisteNum(J,listeRpzOrbiteOTC,VK0t,K_K0tD0);print("listeNum : done");

	print("\n","beginning : ThetaPolynomial");
	listeTFE=vector(#listeFieldArchiOTC,i,0);
	failorbit=vector(#listeFieldArchiOTC,i,i);
	listeThetaPolynomial = listeThetaPolynomialCore(J,listeRpzGx,listeNum,listeOrbiteArchiAlpha,listeFieldArchiOTC,HK0tD0,K_K0tD0,VK0tD0,HK0t,VK0t,correspondance,listeTFE,failorbit);
	\\listeThetaPolynomial=vector(#listeRpzOrbiteOTC,i,ThetaPolynomial(J,listeRpzOrbiteOTC[i],listeRpzGx[i],HK0tD0,VK0t,HK0t,K_K0tD0,listeFieldArchiOTC[i],listeNum[i],listeOrbiteArchiAlpha[i]));
	print("end");
	return(listeThetaPolynomial);
}
*/






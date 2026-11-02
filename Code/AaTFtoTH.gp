\\ \r AaTechnique.gp
\\ \r AaMakdisiMachinery.gp
\\ \r AaTP.gp

\\From Theta Fractions to Theta Hyperplanes

\\********************************************* INTERPOLATION LAMBDA *****************************************************
\\************************************************************************************************************************

\\****************************************** ETC Lambda ******************************************************************

\\WTC representing a ETC divisor beta
\\we check if there is a (small) d such that beta + dP0 >=0 is effective.
\\generically, it should work for d=1



\\****************************************** General Lambda **************************************************************


InitWArfGAdic(J,CCurve,P0,flagETCeffective)=
\\TODO we may avoid to compute WArfETCnonVanish in all the cases where there is no nonèvanishing ETC. But it would require to add some flags
{
	my(WArfGAdic,WArfETC,WArfETCnonVanish,dim0);
	dim0=getd0(J) +getg(J) -1;
	WArfGAdic=PoltoMakdisi(J,RiemannRochTM(CCurve,[P0,dim0]));
	if(flagETCeffective,
		WArfETC=PoltoMakdisi(J,RiemannRochTM(CCurve,[P0,dim0+1]));
		WArfETCnonVanish=PoltoMakdisi(J,RiemannRochTM(CCurve,[P0,dim0-1]));
		WArfGAdic=[WArfGAdic,WArfETC,WArfETCnonVanish];
		);
	return(WArfGAdic);
}

WArfGAdicflagAlpha(J,P0,WArfGAdicInit,flagETCeffective,flagETCnonVanish)=
\\we always have WArfAlpha=WArfGAdic[4] when we compute it
{
	my(WArfAlpha,WArfGAdic,WArfETCnonVanish,dim0);
	dim0=getd0(J) +getg(J) -1;
	WArfGAdic=WArfGAdicInit;
	if((flagETCeffective&&flagETCnonVanish)==0,  \\a OTC with a flagAlpha=0 has to be dealt with WArfETCnonVanish (position 3)
		WArfETCnonVanish=PoltoMakdisi(J,RiemannRochTM(CCurve,[P0,dim0-1]));
                if(type(WArfGAdic)=="t_Mat",
			WArfGAdic=[WArfGAdic,0,WArfETCnonVanish];
			);
		if((type(WArfGAdic)=="t_VEC")&&(#WArfGAdic==2),
			WArfGAdic=concat(WArfGAdic,[WArfETCnonVanish]));
		);

	if(flagETCeffective&&flagETCnonVanish,
		print("check RR dim pour WArfAlpha: d0+g-3=",dim0-2);
		WArfAlpha=PoltoMakdisi(J,RiemannRochTM(CCurve,[P0,dim0-2]));\\for ETC vanishing with flagAlpha=0
		WArfGAdic=WArfGAdicInit;
		if(type(WArfGAdic)=="t_MAT",WArfGAdic=[WArfGAdic,0,0,WArfAlpha]);
		if((type(WArfGAdic)=="t_VEC")&&(#WArfGAdic==2),WArfGAdic=concat(WArfGAdic,[0,WArfAlpha]));
		if((type(WArfGAdic)=="t_VEC")&&(#WArfGAdic==3),WArfGAdic=concat(WArfGAdic,[WArfAlpha]));
		);
	
	return(WArfGAdic);
}


CheckBetaEffective(J,WTC,HDp,WArfGAdic,dimarff,parity,vanishing)=
{
	my(d0,g,dp,JFF,WFF,HDpFF,WArf,HDFF,dim);
	d0=getd0(J);
	g=getg(J);
	dp=d0-g;
	JFF=RedJ(J,1);
	WFF=RedPointJ(J,WTC,1);
	HDpFF=RedPointJ(J,HDp,1);
	WArf=RedPointJ(J,WArfGAdic,1);
	if(default(debug)>2,print("dim Arf = ",#WArf));	
	
	dim=3*d0-(dp+1)-g+1;
	if(parity==0,
        	if(vanishing==0,dim=dim+1,
                	dim=dim-1);
		);
	\\we have HDFF=H°(B) with B~2D0-dP0+D, even when flagAlpha=0
	HDFF=DivAddTM(JFF,HDpFF,WArf,dim);\\H°(4D0-dP0-E1-E) for OTC
	if(default(debug)>2,print("actual dim = ",#HDFF," = ",dim));
        HDFF=DivSubDimKFF(JFF,HDFF,getWn(JFF,2),getKn(JFF,2));\\H°( (2D0-E1-E) -dP0 ) for OTC      ~ D-dP0 ~ Beta
	if(HDFF==0,print("\n error DivSub CheckBetaEffective \n"));
	return(#HDFF);
}

/*
InitdimArf(J)=
\\give the dimension to set in the divAdd of CheckBetaEffective if there is no
{
	my(dimArf,d0,g,dp,dim0,dim);
	dimArf=[0,0,0,0];
	d0=getd0(J);
	g=getg(J);
        dp=d0-g;
	dim0=3*d0-(dp+1)-g+1;\\classic case
	dim=dim0;
	dimArf[1]=dim;
	dim=dim0-1;
	dimArf[2]=dim;
	dim=dim0+1;
	dimArf[3]=dim;
	dim=dim0-2;
	dimArf[4]=dim;
	return(dimArf);
}
*/

InitLambdaBetaGAdic(J,WTC,WArfGAdic,dimarff,parity,vanishing,Gx,flagAlpha,W2plusP0)=
\\parity is 1 or 0 according to the parity of WTC
\\WArf=H°(2D0 -dP0) = H°(D0 +(g-1)P0) for OTC		position 1
\\WArf=H°(2D0 -dP0 +P0) = H°(D0 +gP0) for ETC (vanishing=1)		position 2
\\WArf=H°(2D0 -dP0 -P0) = H°(D0 +(g-2)P0) for ETC (vanishing=0)		position 3
\\if non classic Alpha, we take WArf with -P0 than the usual one
\\WArf=H°(2D0 -dP0 -2P0)		for ETC (vanishing=0 flagAlpha=0)	position 4
\\d=d0 -g +1
\\position corresponds to the position in the global WArfGAdic list
{
	my(V2,WE1,d0,g,dp,d,Wsp,Ws,s,sVp,HD,LambdaBetaGAdic,HDp,dimBeta,flag,Warf,dim,K);
	V2=getWn(J,2);
	WE1=getWE(J)[1][1];\\J[14][1][1];
	d0=getd0(J);
	g=getg(J);
	dp=d0-g;
	flag=1;

	\\print("1");
/*
\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\We are recomputing listeRpzGx
	Wsp=DivAddTM(J,WTC,WE1,d0+(2*d0-dp)-g+1); \\H°(4D0-E1-D)
	print("long? it may be necessary to try with another WE1/WE2 (P1 and/or P2 to change)");
	Ws=DivSubTM(J,Wsp,V2,1,getKn(J,2)); \\H°(2D0-E1-D)
	if(Ws==0,
		print("\n","WE1 and/or WE2 fail. Change P1 and/or P2","\n");
		flag=0;
		);
	s=Ws[,1];\\(s)=D+E1+E-2D0
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*/
	s=Gx;
	
	\\print("3");
	if(flagAlpha==0,V2=W2plusP0);
	sVp=mulPointWise(J,s,V2);\\H°(4D0-D-E1-E)
	if(default(debug)>2,print("dim sVp = ",#sVp," = ",2*d0-g+1-(flagAlpha-1)));
	\\sVp=mulPointWise(J,s,WArf);\\H°(4D0-dP0-D-E1-E)
	dim=d0-g+1;
	K=getKn(J,2);
	if(flagAlpha==0,
		dim=dim+1;
		K=getKn(J,3);
		);
	HDp=DivSubTM(J,sVp,WTC,dim,K);\\H°(2D0-E1-E)
	if(HDp==0,print("\n error DivSub InitLambdaBetaGAdic \n"));

	\\print("4");
	dimBeta=CheckBetaEffective(J,WTC,HDp,WArfGAdic,dimarff,parity,vanishing);
		\\WTC representing a ETC divisor beta
		\\we check if there is a (small) d such that beta + dP0 >=0 is effective.
		\\generically, it should work for d=1
	if(default(debug)>2,print("parity: ",parity," Beta (+/-dP0) uniquely Effective? ",dimBeta==1));
	LambdaBetaGAdic=1;
	\\print("5");
	if(dimBeta==1,
		dim=3*d0-(dp+1)-g+1;
		if(parity==0,
			if(vanishing==0,dim=dim+1,
				dim=dim-1);
			);
		HD=DivAddTM(J,HDp,WArfGAdic,dim);\\H°(4D0-dP0-E1-E)  for OTC
		HD=DivSubTM(J,HD,V2,1,getKn(J,2));\\H°( (2D0-E1-E) -dP0 )   for OTC
		\\HD=DivSubTM(J,sVp,WTC,1,getKn(J,2));\\H°(2D0-dP0-E1-E)=H°(Dbeta'-dP0)
		if(HD!=0,s=HD[,1],
			print("\n Error DivSub InitLambdaBetaGAdic \n");
			s=0;
			);
		LambdaBetaGAdic=s;
		);
	\\if(LambdaBetaGAdic==0,print("\n","error LambdaBetaGAdic","\n"));
	return([dimBeta,LambdaBetaGAdic,flag]);
}


InitlisteOrbiteLambdaGAdic(J,listeFieldArchiOTC,correspondance,listeRpzLambdaGAdic)=
\\GAdic rpz of the Global Q-Galois orbits
{
	my(listeOrbiteLambdaGAdic,lindex,OrbiteLambda,flagLambda,Lambda);
	listeOrbiteLambdaGAdic=vector(#listeFieldArchiOTC,i,0);
	for(i=1,#listeFieldArchiOTC,
                lindex=correspondance[i];\\the indexes of the p-adic orbits corresponfing to the Galois orbit indexed by i
		\\Gathering the Lambda in each p-adic orbit of the same Galois orbit
		OrbiteLambda=vector(#lindex,k,0);
                flagLambda=1;
                for(k=1,#lindex,
                        if(flagLambda,
                                [flagLambda,Lambda]=listeRpzLambdaGAdic[lindex[k]];
                                if(flagLambda,
                                        OrbiteLambda[k]=Lambda);
                                );
                        );
		listeOrbiteLambdaGAdic[i]=[flagLambda,OrbiteLambda];
		);
	return(listeOrbiteLambdaGAdic);
}
/*
	my(listeLambdaRpzGAdic,lindex,WTC,dimBeta,flagLambda,OrbiteLambdaGAdic);
        listeLambdaRpzGAdic=vector(#listeFieldArchiOTC,i,0);
        for(i=1,#listeFieldArchiOTC,
                lindex=correspondance[i];\\the indexes of the p-adic orbits corresponfing to the Galois orbit indexed by i

                \\Gathering the Lambda in each p-adic orbit of the same Galois orbit
                WTC=listeRpzOrbiteOTC[lindex[1]];
                print("Lamba: new Global orbit (",i,"/",#listeFieldArchiOTC,")");
                [dimBeta,LambdaBetaRpzGAdic]=InitLambdaBetaGAdic(J,WTC,WArfGAdic);
                print("dim Beta = ",dimBeta);
                if(dimBeta==1,
                        flagLambda=1;
                        OrbiteLambdaGAdic=vector(#lindex,i,0);
                        OrbiteLambdaGAdic[1]=LambdaBetaRpzGAdic;
                        for(k=2,#lindex,
                                OrbiteLambdaGAdic[k]=InitLambdaBetaGAdic(J,listeRpzOrbiteOTC[lindex[k]],WArfGAdic)[2];
                                );
                        );

                if(dimBeta!=1,
                        flagLambda=0;
                        OrbiteLambdaGAdic=0;
                        );

                listeLambdaRpzGAdic[i]=[flagLambda,OrbiteLambdaGAdic];
                );
        return(listeLambdaRpzGAdic);
*/


InitlisteLambdaArchi(J,listeOrbiteLambdaGAdic,correspondance,listeOrbiteArchiAlpha,listeFieldArchiOTC,failorbit)=
{
        my(W2,H2D0,flagLambda,listeLambdaArchi,lindex,OrbiteArchiAlpha,Field,OrbiteLambdaGAdic,OrbiteLambdaCoefBasis,LambdaArchi,flagInterpolation);
        W2=getWn(J,2);
        \\H2D0=J[4][2];
	H2D0=getW2Pol(J);
	\\listeLambdaRpzGAdic=InitlisteLambdaRpzGAdic(J,listeFieldArchiOTC,correspondance,listeRpzOrbiteOTC,WArfGAdic);
        listeLambdaArchi=vector(#listeFieldArchiOTC,i,0);
        for(i=1,#listeFieldArchiOTC,
		if(IsInList(failorbit,i),
                	lindex=correspondance[i];\\the indexes of the p-adic orbits corresponfing to the Galois orbit indexed by i
                	OrbiteArchiAlpha=listeOrbiteArchiAlpha[i];
                	Field=listeFieldArchiOTC[i];
			[flagLambda,OrbiteLambdaGAdic]=listeOrbiteLambdaGAdic[i];
                	if(flagLambda==1,
                        	OrbiteLambdaCoefBasis=vector(#lindex,k,ColtoPolCoefClean(J,OrbiteLambdaGAdic[k],W2,0)[2]);
                        	[flagInterpolation,LambdaArchi]=InterpolationTheta(J,OrbiteLambdaCoefBasis,Field,H2D0,OrbiteArchiAlpha);
                        	if(flagInterpolation==0,
                                	LambdaArchi=0;
                                	print("InitlisteLambdaArchi ERROR: InterpolationTheta failed");
                                	);
                        	);
                	if(flagLambda!=1,
                        	LambdaArchi=1;
                        	);
                	listeLambdaArchi[i]=[flagLambda,LambdaArchi];
                	);
		);
        return(listeLambdaArchi);
}

/*
InitlisteLambdaArchi(J,correspondance,listeRpzOrbiteOTC,listeOrbiteArchiAlpha,listeFieldArchiOTC,WArfGAdic)=
{
	my(W2,H2D0,flagLambda,listeLambdaArchi,lindex,OrbiteArchiAlpha,Field,WTC,dimBeta,LambdaBetaRpzGAdic,OrbiteLambdaGAdic,OrbiteLambdaCoefBasis,LambdaArchi,flagInterpolation);
	W2=getWn(J,2);
	H2D0=J[4][2];
	listeLambdaArchi=vector(#listeFieldArchiOTC,i,0);
	for(i=1,#listeFieldArchiOTC,
		lindex=correspondance[i];\\the indexes of the p-adic orbits corresponfing to the Galois orbit indexed by i
		OrbiteArchiAlpha=listeOrbiteArchiAlpha[i];
                Field=listeFieldArchiOTC[i];

		\\Gathering the Lambda in each p-adic orbit of the same Galois orbit
		WTC=listeRpzOrbiteOTC[lindex[1]];
		print("Lamba: new Global orbit (",i,"/",#listeFieldArchiOTC,")");
		[dimBeta,LambdaBetaRpzGAdic]=InitLambdaBetaGAdic(J,WTC,WArfGAdic);
		print("dim Beta = ",dimBeta);
		if(dimBeta==1,
			flagLambda=1;
			OrbiteLambdaGAdic=vector(#lindex,i,0);
			OrbiteLambdaGAdic[1]=LambdaBetaRpzGAdic;
			for(k=2,#lindex,
				OrbiteLambdaGAdic[k]=InitLambdaBetaGAdic(J,listeRpzOrbiteOTC[lindex[k]],WArfGAdic)[2];
				);

			OrbiteLambdaCoefBasis=vector(#lindex,k,ColtoPolCoefClean(J,OrbiteLambdaGAdic[k],W2));
			[flagInterpolation,LambdaArchi]=InterpolationTheta(J,OrbiteLambdaCoefBasis,Field,H2D0,OrbiteArchiAlpha);
			if(flagInterpolation==0,
				LambdaArchi=0;
				print("InitlisteLambdaArchi ERROR: InterpolationTheta failed");
				);
			);

		if(dimBeta!=1,
			flagLambda=0;
			LambdaArchi=1;
			);

		listeLambdaArchi[i]=[flagLambda,LambdaArchi];
		);
	return(listeLambdaArchi);
}
*/

listeLambdaArchiMain(J,listeFieldArchiOTC,correspondance,listeRpzLambdaGAdic,listeOrbiteArchiAlpha,failorbit,flagTC,flagETCeffective,parity)=
{
	my(listeOrbiteLambdaGAdic,listeLambdaArchi,i,lLambda);
	if(flagTC==0,
		listeOrbiteLambdaGAdic=InitlisteOrbiteLambdaGAdic(J,listeFieldArchiOTC,correspondance,listeRpzLambdaGAdic);
        	listeLambdaArchi=InitlisteLambdaArchi(J,listeOrbiteLambdaGAdic,correspondance,listeOrbiteArchiAlpha,listeFieldArchiOTC,failorbit);
		);
	if(flagTC,
		if(parity==1,i=1);
		if(parity==0,i=2);
		if(flafETCeffective==0,lLambda=listeRpzLambdaGAdic,
			lLambda=listeRpzLambdaGAdic[i]);
		listeOrbiteLambdaGAdic=InitlisteOrbiteLambdaGAdic(J,listeFieldArchiOTC[i],correspondance[i],lLambda);
        	listeLambdaArchi=InitlisteLambdaArchi(J,listeOrbiteLambdaGAdic,correspondance[i],listeOrbiteArchiAlpha[i],listeFieldArchiOTC[i],failorbit[i]);
		);
	return(listeLambdaArchi);
}



MakeTFEffective(listeThetaPolynomialInit,failorbit,listeLambdaArchi,PolCurve)=
\\EX InitlisteTFEffective(PolCurve,listeThetaPolynomialInit,listeLambdaArchi,failorbit)=
\\Since listeThetaPolynomialInit comes from Effective's technique, if it is not empty, it is only TFE
\\Becareful: indexes not in failorbit are supposed to refer to TF that are effective 
{
	my(listeIndexTFE,listeThetaFraction,flagLambda,flagTheta);
	listeIndexTFE=[];
	listeThetaFraction=listeThetaPolynomialInit;
	for(i=1,#listeThetaFraction,
		if(IsInList(failorbit,i)==0,listeIndexTFE=concat(listeIndexTFE,[i]),
			flagLambda=listeLambdaArchi[i][1];
			if(flagLambda==0,listeTFE[i]=[0,0],
				flagTheta=listeThetaFraction[i][1];
				if(flagTheta,
					print("\n","new TF made Effective in TWO PHASES");
					TFE=lift(Mod(listeThetaFraction[i][2]*listeLambdaArchi[i][2]^2,PolCurve));
					listeThetaFraction[i][2]=TFE;
					listeIndexTFE=concat(listeIndexTFE,[i]);
					);
				);
			);
		);
	return([listeThetaFraction,listeIndexTFE]);
}

MakeTFEffectiveMain(J,PolCurve,listeThetaPolynomialInit,listeFieldArchiOTC,correspondance,listeRpzLambdaGAdic,listeOrbiteArchiAlpha,failorbit,flagTC,flagETCeffective)=
\\Old InitlisteTFEffective(PolCurve,listeThetaPolynomialInit,listeLambdaArchi,failorbit)=
\\Since listeThetaPolynomialInit comes from Effective's technique, if it is not empty, it is only TFE
\\Becareful: indexes not in failorbit are supposed to refer to TF that are effective 
{
	my(listeLambdaArchiOdd,listeThetaFractionOdd,listeIndexTFEOdd,listeThetaFraction,listeIndexTFE,listeLambdaArchiEven,listeThetaFractionEven,listeIndexTFEEven);

	\\************************* OTC *************************
	listeLambdaArchiOdd=listeLambdaArchiMain(J,listeFieldArchiOTC,correspondance,listeRpzLambdaGAdic,listeOrbiteArchiAlpha,failorbit,flagTC,flagETCeffective,1);
	if(flagTC,[listeThetaFractionOdd,listeIndexTFEOdd]=MakeTFEffective(listeThetaPolynomialInit[1],failorbit[1],listeLambdaArchiOdd,PolCurve),
		[listeThetaFraction,listeIndexTFE]=MakeTFEffective(listeThetaPolynomialInit,failorbit,listeLambdaArchiOdd,PolCurve));

	if(flagTC,
		\\***************** ETC ************************
		if(flagETCeffective==0,
			listeThetaFractionEven=listeThetaPolynomialInit[2];
			listeIndexTFEEven=0;
			);
		if(flagETCeffective,
			listeLambdaArchiEven=listeLambdaArchiMain(J,listeFieldArchiOTC,correspondance,listeRpzLambdaGAdic,listeOrbiteArchiAlpha,failorbit,flagTC,flagETCeffective,0);
			[listeThetaFractionEven,listeIndexTFEEven]=MakeTFEffective(listeThetaPolynomialInit[2],failorbit[2],listeLambdaArchiEven,PolCurve);
			);	
		listeThetaFraction=[listeThetaFractionOdd,listeThetaFractionEven];
		listeIndexTFE=[listeIndexTFEOdd,listeIndexTFEEven];
		);
	print("listeIndexTFE = ",listeIndexTFE);
	return([listeThetaFraction,listeIndexTFE]);
}





\\************************************************************ TFE in one phase ***************************************************************
/*
InversesThetaInvMakdisi(J,sThetaInv)=
{
	my(wrongindex,sTheta);
	wrongindex=[];
	sTheta=vector(#sThetaInv,k,0)~;
	for(k=1,#sTheta,
		\\print("k= ",k);
		\\we test if sThetaInv[k] is inversible. So we look if its reduction is 0. Since p is of good reduction, we check p/lift(a)
		if(Mod(sThetaInv[k],getpMak(J))!=0,sTheta[k]=ZpXQ_inv(sThetaInv[k],getTNicolas(J),getpMak(J),geteJ(J)),	\\sTheta[k]=divpadicElt(J,1,sThetaInv[k]),
			wrongindex=concat(wrongindex,[k]);
			sTheta[k]=+oo;
			);
		);
	return([sTheta,wrongindex]);
}

InitlOrbsTheta(J,lOrbsThetaInv)=
{
	my(lOrbsTheta,wrongindex,sTheta,wi);
	lOrbsTheta=vector(#lOrbsThetaInv,k,0);
        wrongindex=[];
        for(k=1,#lOrbsTheta,
        	[sTheta,wi]=InversesThetaInvMakdisi(J,lOrbsThetaInv[k]);
                wrongindex=concat(wrongindex,wi);
                lOrbsTheta[k]=sTheta;
                );
	return([lOrbsTheta,wrongindex]);
}

InitOrbiteTFEGAdic(J,wrongindex,OrbiteLambdaGAdic,lOrbsTheta)=
{
	my(OrbiteTFEGAdic,sTheta,lambda,c);
	OrbiteTFEGAdic=vector(#OrbiteLambdaGAdic,k,0);
	for(k=1,#OrbiteTFEGAdic,
		sTheta=lOrbsTheta[k];
		lambda=OrbiteLambdaGAdic[k];
		c=vector(#sTheta,l,0)~;
		for(l=1,#c,
			if(IsInList(wrongindex,l),c[l]=+oo,
				\\print("sTheta == 1 ? ",sTheta[l]==sTheta[1]);
				c[l]=mulpadicElt(J,sTheta[l],mulpadicElt(J,lambda[l],lambda[l])); \\f*lambda^2
				\\print("lambda nul ? ",lambda[l]==0);
				\\print("c == 1 ? ",c[l]==c[1]);
				);
			);
		OrbiteTFEGAdic[k]=c;
		);
	return(OrbiteTFEGAdic);
}

InitOrbiteTFECoefBasis(J,OrbiteTFEGAdic,WK0,wrongindex)=
{
	my(OrbiteTFECoefBasis,flagwi);
	OrbiteTFECoefBasis=vector(#OrbiteTFEGAdic,k,0);
	flagwi=1;
	for(k=1,#OrbiteTFEGAdic,
		if(flagwi,
			[flagwi,CoefBasis]=ColtoPolCoefClean(J,OrbiteTFEGAdic[k],WK0,wrongindex);	
			);
		OrbiteTFECoefBasis[k]=CoefBasis;
		);
	return([flagwi,OrbiteTFECoefBasis]);
}

InitlisteRpzTFECoefBasis(J,listeRpzLambdaGAdic,listeRpzGx,listeNum,WK0)=
{
	my(listeRpzTFECoefBasis,flagLambda,LambdaGAdic,flagTFE,TFECoefBasis,ThetaInv,sTheta,wrongindex,TFEGAdic,flagwi);
	listeRpzTFECoefBasis=vector(#listeLambdaRpzGAdic,i,0);
	for(i=1,#listeRpzLambdaGAdic,
		[flagLambda,LambdaGAdic]=listeRpzLambdaGAdic[i];
                flagTFE=1;      
                if(flagLambda==0,
                        print("dim Beta != 1. orbit #",i);
                        flagTFE=0;
			TFECoefBasis=[];
			);
                if(flagLambda==1,
                        ThetaInv=InitsThetaDivSubInv(J,listeRpzGx[i],listeNum[i]);
                        [sTheta,wrongindex]=InitTheta(J,ThetaInv);
                        \\OrbiteTFEGAdic=InitOrbiteTFEGAdic(J,wrongindex,OrbiteLambdaGAdic,lOrbsTheta);\\print("c == 1 ? ",OrbiteTFEGAdic[1][1]==OrbiteTFEGAdic[1][2]);
                        TFEGAdic=InitTFEGAdic(J,wrongindex,LambdaGAdic,sTheta);
			\\[flagwi,OrbiteTFECoefBasis]=InitOrbiteTFECoefBasis(J,OrbiteTFEGAdic,WK0,wrongindex);\\print("CoefBasis = ",OrbiteTFECoefBasis[1]);
                        [flagwi,TFECoefBasis]=InitTFECoefBasis(J,TFEGAdic,WK0,wrongindex);
			if(flagwi==0,
                                print("we couldn't identify the TFE because too much poles on the Makdisi evaluation Points Z. orbit #",i);
                                flagTFE=0;
	P0p=[0,1,1];
				TFECoefBasis=[];
				);
			);
		listeRpzTFECoefBasis[i]=[flagTFE,TFECoefBasis];
		);
	return(listeRpzTFECoefBasis);
}
*/

InitlisteTFE(J,listeFieldArchiOTC,listeRpzTFECoefBasis,listeOrbiteArchiAlpha,correspondance,HK0)=
\\We try to interpolate the TFE right away as such in only one phase
{
	my(failorbit,listeTFE,lindex,OrbiteArchiAlpha,Field,flagTFE,OrbiteTFECoefBasis,flagCoefBasis,TFECoefBasis,flagInterpolation,TFE);
	if(default(debug)>1,print("we try to interpolate the TFE in one phase"));
        \\listeLambdaRpzGAdic=InitlisteLambdaRpzGAdic(J,listeFieldArchiOTC,correspondance,listeRpzOrbiteOTC,WArfGAdic);
        listeTFE=vector(#listeFieldArchiOTC,i,0);
	failorbit=[];
        for(i=1,#listeFieldArchiOTC,
                lindex=correspondance[i];\\the indexes of the p-adic orbits corresponding to the Galois orbit indexed by i
                OrbiteArchiAlpha=listeOrbiteArchiAlpha[i];
                Field=listeFieldArchiOTC[i];
		
		OrbiteTFECoefBasis=vector(#lindex,k,0);
		flagCoefBasis=1;
		for(k=1,#lindex,
			if(flagCoefBasis,
				[flagCoefBasis,TFECoefBasis]=listeRpzTFECoefBasis[lindex[k]];
				if(flagCoefBasis,
					OrbiteTFECoefBasis[k]=TFECoefBasis);
				);
			);
		flagTFE=flagCoefBasis;

		\\[flagTFE,OrbiteTFECoefBasis]=listeOrbiteTFECoefBasis[
		/*
                [flagLambda,OrbiteLambdaGAdic]=listeLambdaRpzGAdic[i];
		flagTFE=1;
                if(flagLambda==0,
			print("dim Beta != 1. orbit #",i);
			flagTFE=0;);
		if(flagLambda==1,
			lOrbsThetaInv=vector(#lindex,k,InitsThetaDivSubInv(J,listeRpzGx[lindex[k]],listeNum[lindex[k]]));
			[lOrbsTheta,wrongindex]=InitlOrbsTheta(J,lOrbsThetaInv);\\print("sTheta == 1 ? ",lOrbsTheta[1][1]==lOrbsTheta[1][2]);
			OrbiteTFEGAdic=InitOrbiteTFEGAdic(J,wrongindex,OrbiteLambdaGAdic,lOrbsTheta);\\print("c == 1 ? ",OrbiteTFEGAdic[1][1]==OrbiteTFEGAdic[1][2]);
			[flagwi,OrbiteTFECoefBasis]=InitOrbiteTFECoefBasis(J,OrbiteTFEGAdic,WK0,wrongindex);\\print("CoefBasis = ",OrbiteTFECoefBasis[1]);
			if(flagwi==0,
				print("we couldn't identify the TFE because too much poles on the Makdisi evaluation Points Z. orbit #",i);
				flagTFE=0;);
*/
		if(flagTFE,
			[flagInterpolation,TFE]=InterpolationTheta(J,OrbiteTFECoefBasis,Field,HK0,OrbiteArchiAlpha);
			\\print("flag = ",flagInterpolation," TFE = ",TFE);
                	if(flagInterpolation==0,
				print("Interpolation failed. orbit #",i);
				flagTFE=0;);
                        );
                if(flagTFE==0,
                        TFE=1;
			failorbit=concat(failorbit,[i]);
                        );
                listeTFE[i]=[flagTFE,TFE];
                );
	print("         failorbit: ",failorbit);
        return([listeTFE,failorbit]);
}










\\************************************************************ THETA HYPERPLANES *************************************************************

FromTFtoTH(TF,Omega,Field)=
{
	my(zF,bas,dimZ,deg,dimc,Z,M,f,Y,c,cpol,TH);
	zF=variable(Field.pol);
	bas=Omega;
	dimZ=50*#bas;
	deg=poldegree(Field.pol);
	dimc=#bas*deg;
	Z=vector(dimZ,i,[random(100*dimc),random(100*dimc),random(100*dimc)]);
	M=matrix(dimZ,dimc,i,j,0);
	for(l=1,dimZ,
		for(i=1,#bas,
			for(j=1,deg,
				M[l,(i-1)*deg+j]=substvec(z^(j-1)*bas[i],[x,y,z],Z[l]);
				);););
	f=lift(TF);
	f=subst(f,zF,z);
	Y=vector(dimZ,l,substvec(f,[x,y,z],Z[l]))~;
	c=matinverseimage(M,Y);
	TH=vector(#bas,i,0);
	for(i=1,#bas,
		cpol=0;
		for(j=1,deg,
			cpol=cpol+ c[(i-1)*deg+j]*z^(j-1);
			);
		TH[i]=Mod(cpol,subst(Field.pol,zF,z));
		);
	return(TH);
}


InitlisteTH(listeThetaFraction,Omega,listeFieldArchiOTC,listeIndexTFE)=
{
	my(listeTH);
	listeTH=vector(#listeThetaFraction,i,0);
	for(i=1,#listeTH,
		if(IsInList(listeIndexTFE,i),
			listeTH[i]=FromTFtoTH(listeThetaFraction[i][2],Omega,listeFieldArchiOTC[i]);
			);
		);
	return(listeTH);
}

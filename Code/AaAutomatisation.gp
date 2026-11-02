\\******************************************************************************************
\\*                                     Automatisation                                     *
\\******************************************************************************************





\\***************************************** Chi ************************************************************
\\**********************************************************************************************************

\\flagmethod = "Nicolas", "Classic", "TF", "TFwithCL"
\\flagTC=1 if we want to compute all the TC. flagTC=0 we compute just the OTC

ComputeChiOTCCore(J,listeRpzOrbiteOTC,listeOrbiteOdd,flagmethod,flagTC,CCurve,P0,P1,P2,Frob)=
{
	my(Chi,lOrbiteAlpha,flagChi);
	if(flagmethod=="Classic",
                if(flagTC,[Chi,lOrbiteAlpha,flagChi]=InitClassicChiOTC(J,listeRpzOrbiteOTC[1],listeOrbiteOdd[1]),
                        [Chi,lOrbiteAlpha,flagChi]=InitClassicChiOTC(J,listeRpzOrbiteOTC,listeOrbiteOdd)));

        if(flagmethod=="Nicolas",
                if(flagTC,[Chi,lOrbiteAlpha,flagChi]=InitChiOTCNicolas(J,listeRpzOrbiteOTC[1],listeOrbiteOdd[1],CCurve,P0,P1,P2,Frob),
                        [Chi,lOrbiteAlpha,flagChi]=InitChiOTCNicolas(J,listeRpzOrbiteOTC,listeOrbiteOdd,CCurve,P0,P1,P2,Frob)));
	return([Chi,lOrbiteAlpha,flagChi]);
}

ComputeChiOTC(J,listeRpzOrbiteOTC,listeOrbiteOdd,flagmethod,flagTC,CCurve,P0,P1,P2,Frob)=
{
	my(Chi,lOrbiteAlpha,flagChi);
	\\J=PrepChi(J,CCurve,P0,P1,P2,Frob);
        [Chi,lOrbiteAlpha,flagChi]=ComputeChiOTCCore(J,listeRpzOrbiteOTC,listeOrbiteOdd,flagmethod,flagTC,CCurve,P0,P1,P2,Frob);
	/*
	if((flagChi==0)&&(P2!=0)&&(P1!=0),
		print("trying another WE1/WE2");
                J=PrepChi(J,CCurve,P0,P1,0,Frob);
                [Chi,lOrbiteAlpha,flagChi]=ComputeChiOTCCore(J,listeRpzOrbiteOTC,listeOrbiteOdd,flagmethod,flagTC,CCurve,P0,P1,P2,Frob);
		);
        if((flagChi==0)&&(P2!=0),
		print("trying another WE1/WE2");
                J=PrepChi(J,CCurve,P0,P2,0,Frob);
                [Chi,lOrbiteAlpha,flagChi]=ComputeChiOTCCore(J,listeRpzOrbiteOTC,listeOrbiteOdd,flagmethod,flagTC,CCurve,P0,P1,P2,Frob);
		);
        if((flagChi==0)&&(P1!=0),
		print("trying another WE1/WE2");
                J=PrepChi(J,CCurve,P0,0,0,Frob);
                [Chi,lOrbiteAlpha,flagChi]=ComputeChiOTCCore(J,listeRpzOrbiteOTC,listeOrbiteOdd,flagmethod,flagTC,CCurve,P0,P1,P2,Frob);
		);
	*/
	return([J,Chi,lOrbiteAlpha,flagChi]);
}

UpdateChiOTC(ChiOTCInit,ChiETC,listeOrbiteAlphaInit,listeOrbiteAlphaEven,flagChiOTCInit,flagChiETC)=
{
	my(ChiOTC,listeOrbiteAlpha,flagChiOTC);
	ChiOTC=[ChiOTCInit,ChiETC];
        listeOrbiteAlpha=[listeOrbiteAlphaInit,listeOrbiteAlphaEven];
        flagChiOTC=[flagChiOTCInit,flagChiETC];
	return([ChiOTC,listeOrbiteAlpha,flagChiOTC]);
}

InitChiETCAlpha(JInit,listeRpzOrbiteOTC,listeOrbiteOdd,flagmethodTC,CCurve,P0,P1,P2,Frob)=
{
	my(J,ChiETC,listeOrbiteAlphaEven,flagChiETC);
	J=JInit;
	if((flagmethodTC=="Nicolas")||(flagmethodTC=="Classic"),
                        print("         Computation ChiETC by alpha methods");
                        [J,ChiETC,listeOrbiteAlphaEven,flagChiETC]=ComputeChiOTC(J,listeRpzOrbiteOTC[2],listeOrbiteOdd[2],flagmethodTC,0,CCurve,P0,P1,P2,Frob);
                        );
	\\ChiOTC=[ChiOTCInit,ChiETC];
        \\listeOrbiteAlpha=[listeOrbiteAlphaInit,listeOrbiteEvenAlpha];
        \\flagChiOTC=[flagChiOTCInit,flagChiETC];
	return([J,ChiETC,listeOrbiteAlphaEven,flagChiETC]);
}

InitChiOTCAlphaMain(JInit,listeRpzOrbiteOTC,listeOrbiteOdd,flagmethod,flagTC,CCurve,P0,P1,P2,Frob)=
\\according to flagmethod, listeRpzOrbiteOTC or listeRpzTFECoefBasis can be set to 0
{
	my(J,ChiOTC,listeOrbiteAlpha,flagChiOTC,ChiETC,listeOrbiteEvenAlpha,flagChiETC);
	J=JInit;
	ChiOTC=0;
	listeOrbiteAlpha=0;
	flagChiOTC=0;
	[J,Chi,lOrbiteAlpha,flagChi]=ComputeChiOTC(J,listeRpzOrbiteOTC,listeOrbiteOdd,flagmethod,flagTC,CCurve,P0,P1,P2,Frob);
	\\if(flagTC,[J,ChiOTC,listeOrbiteAlpha,flagChiOTC]=ComputeChiOTC(J,listeRpzOrbiteOTC[1],listeOrbiteOdd[1],flagmethod,flagTC,CCurve,P0,P1,P2,Frob),
		\\[J,ChiOTC,listeOrbiteAlpha,flagChiOTC]=ComputeChiOTC(J,listeRpzOrbiteOTC,listeOrbiteOdd,flagmethod,flagTC,CCurve,P0,P1,P2,Frob));
	return([J,ChiOTC,listeOrbiteAlpha,flagChiOTC]);
}

InitChiOTCbyTFEMain(J,listeRpzTFECoefBasis,listeOrbiteOdd,flagmethod,flagTC,flagETCeffective,parity)=
{
	my(lCoef,ChiOTC,listeOrbiteAlpha,flagChiOTC);
	if(flagETCeffective==0,lCoef=listeRpzTFECoefBasis,
		lCoef=listeRpzTFECoefBasis[1]);

	if(flagmethod=="TF",
		if(flagTC,[ChiOTC,listeOrbiteAlpha,flagChiOTC]=InitChiOTCbyTF(J,lCoef,listeOrbiteOdd[1]),
			[ChiOTC,listeOrbiteAlpha,flagChiOTC]=InitChiOTCbyTF(J,lCoef,listeOrbiteOdd));
		);

	if(flagmethod=="TFwithCL",
		if(flagTC,[ChiOTC,listeOrbiteAlpha,flagChiOTC]=InitChiOTCbyTFwithCL(J,lCoef,listeOrbiteOdd[1],parity),
			[ChiOTC,listeOrbiteAlpha,flagChiOTC]=InitChiOTCbyTFwithCL(J,lCoef,listeOrbiteOdd,parity));
		);

	return([ChiOTC,listeOrbiteAlpha,flagChiOTC]);
}


InitChiETCbyTFEMain(J,listeRpzTFECoefBasis,listeOrbiteOdd,flagmethodTC,flagETCeffective)=
{
	my(ChiETC,listeOrbiteAlphaEven,flagChiETC);
	if(flagETCeffective==0,print("\n","ERROR InitChiETCbyTFEMain: flagETCeffective should be 1","\n"));
	[ChiETC,listeOrbiteAlphaEven,flagChiETC]=InitChiOTCbyTFEMain(J,listeRpzTFECoefBasis[2],listeOrbiteOdd[2],flagmethodTC,0,0,0);
	\\[ChiOTC,listeOrbiteAlpha,flagChiOTC]=[ChiOTCInit,listeOrbiteAlphaInit,flagChiOTCInit];
	\\if(flagChiETC,
		\\[ChiOTC,listeOrbiteAlpha,flagChiOTC]=UpdateChiOTC(ChiOTCInit,ChiETC,listeOrbiteAlphaInit,listeOrbiteEvenAlpha,flagChiOTCInit,flagChiETC);
		\\);
	return([ChiETC,listeOrbiteAlphaEven,flagChiETC]);
}

InitChiOTCbyTFEMix(JInit,listeRpzOrbiteOTC,listeRpzTFECoefBasis,listeOrbiteOdd,flagmethodbis,flagTC,CCurve,P0,P1,P2,Frob,flagETCeffective)=
{
	my(J,lCoef,ChiOTC,listeOrbiteAlpha,flagChiOTC);
	J=JInit;
	if(flagETCeffective==0,lCoef=listeRpzTFECoefBasis,
                lCoef=listeRpzTFECoefBasis[1]);

	if((flagmethodbis=="Nicolas")||(flagmethodbis=="Classic"),
		[J,ChiOTC,listeOrbiteAlpha,flagChiOTC]=InitChiOTCAlphaMain(J,listeRpzOrbiteOTC,listeOrbiteOdd,flagmethodbis,flagTC,CCurve,P0,P1,P2,Frob));
	if((flagmethodbis=="TFwithCL")||(flagmethodbis=="TF"),
		[ChiOTC,listeOrbiteAlpha,flagChiOTC]=InitChiOTCbyTFEMain(J,lCoef,listeOrbiteOdd,flagmethod,flagTC,listeRpzOrbiteOTC));
	return([J,ChiOTC,listeOrbiteAlpha,flagChiOTC]);
}


\\*********************************************************** CHI J[2] ****************************************

InitChiJ2MainCore(J,flagmethodChiJ2,PicTorsFrobGendata,Auts,Frob,SBasisJ2GAdic)=
{
	my(ChiJ2,listeRpzOrbiteJ2,listeOrbiteJ2,listeOrbiteAlpha,flagdeg,flagroot,flagChiJ2);
	print("Etale Polynomila Chi for J[2]");	
	if(flagmethodChiJ2=="NicolasFromScratch",
		ChiJ2=pictorsgalrep(J,2)[1];);

	if(flagmethodChiJ2=="Nicolas",
		ChiJ2=InitChiJ2Nicolas(J,PicTorsFrobGendata,Auts)[1]);

	if(flagmethodChiJ2=="Classic",
		[listeRpzOrbiteJ2,listeOrbiteJ2]=PrepChiJ2(J,Frob,SBasisJ2GAdic);
        	[ChiJ2,listeOrbiteAlpha]=InitChiJ2(J,listeRpzOrbiteJ2,listeOrbiteJ2);
		);
	print("end");
	return(ChiJ2);
}

InitChiJ2Main(J,flagmethodChiJ2,PicTorsFrobGendata,Auts,Frob,SBasisJ2GAdic,flagTC,ChiOTC)=
{
	my(nbrat,ChiJ2,flagChiJ2,matfact,linearfact);
	nbrat=NbRational(ChiOTC,flagTC);
	matfact=0;
	if((flagTC&&(nbrat>0))==0,ChiJ2=InitChiJ2MainCore(J,flagmethodChiJ2,PicTorsFrobGendata,Auts,Frob,SBasisJ2GAdic),
		ChiJ2=ChiOTC[1]*ChiOTC[2];
		matfact=factor(ChiJ2);
		linearfact=matfact[1,1];
		ChiJ2=ChiJ2/linearfact;
		);
	
	\\ChiJ2=CleanChi(ChiJ2);
        flagChiJ2=(poldegree(gcd(ChiJ2,ChiJ2'))==0);\\flag=1 if everything's okay
        print("Chi J2: ",ChiJ2);
        if(flagChiJ2==0,print("BAD CHI"),
                print("GOOD CHI");
		matfact=factor(ChiJ2);
                PrintOrbitDecomposition(matfact);
                );
	return([ChiJ2,flagChiJ2]);
}





\\*************************************************************************************************************
\\*************************************************************************************************************








\\****************************************** InitlisteRpzLambdaGAdic ******************************************

InitLambdaGAdicCoreBis(J,listeRpzOrbiteOTC,WArfGAdic,dimArf,CCurve,P0,P1,P2,Frob,flagTC,flagETCeffective,parity,flagETCnonVanish,listeRpzGx,W2plusP0)=
{
	my(listeRpzLambdaGAdic,flag);
	\\J=JInit;
	[listeRpzLambdaGAdic,flag]=InitlisteRpzLambdaGAdic(J,listeRpzOrbiteOTC,WArfGAdic,dimArf,flagTC,flagETCeffective,parity,flagETCnonVanish,listeRpzGx,W2plusP0);
	if(flag==0,print("\n error InitlisteRpzLambdaGAdic \n"));
	/*
	if((flag==0)&&(P2!=0)&&(P1!=0),
		print("\n","trying another WE1/WE2","\n");
		J=PrepChi(J,CCurve,P0,P1,0,Frob);
		[listeRpzLambdaGAdic,flag]=InitlisteRpzLambdaGAdic(J,listeRpzOrbiteOTC,WArfGAdic,flagTC,flagETCeffective,parity,flagETCnonVanish);
		);
	if((flag==0)&&(P2!=0),
                print("\n","trying another WE1/WE2","\n");
		J=PrepChi(J,CCurve,P0,P2,0,Frob); 
                [listeRpzLambdaGAdic,flag]=InitlisteRpzLambdaGAdic(J,listeRpzOrbiteOTC,WArfGAdic,flagTC,flagETCeffective,parity,flagETCnonVanish);
                );
	if((flag==0)&&(P1!=0),
		print("\n","trying another WE1/WE2","\n");
		J=PrepChi(J,CCurve,P0,0,0,Frob); 
                [listeRpzLambdaGAdic,flag]=InitlisteRpzLambdaGAdic(J,listeRpzOrbiteOTC,WArfGAdic,flagTC,flagETCeffective,parity,flagETCnonVanish);
                );
	*/
	return([J,listeRpzLambdaGAdic,flag]);
}

InitLambdaGAdicMain(JInit,listeRpzOrbiteOTC,WArfGAdic,dimArf,CCurve,P0,P1,P2,Frob,flagTC,flagETCeffective,flagETCnonVanish,listeRpzGx,W2plusP0)=
{
	my(J,listeRpzLambdaGAdic,flaglisteLambdaGAdic,resETC);
	
	print("\n","If it crashes, please start again with flagEffective=0 or/and flagETCeffective=0 according to the parity of the TC that makes the bug. Maybe start again with flagTFEonephase=0.","\n"); 

	\\OTC
	[J,listeRpzLambdaGAdic,flaglisteLambdaGAdic]=InitLambdaGAdicCoreBis(JInit,listeRpzOrbiteOTC,WArfGAdic,dimArf,CCurve,P0,P1,P2,Frob,flagTC,flagETCeffective,1,flagETCnonVanish,listeRpzGx,W2plusP0);

	if(flagETCeffective,
		\\ETC
		resETC=InitLambdaGAdicCoreBis(J,listeRpzOrbiteOTC,WArfGAdic,dimArf,CCurve,P0,P1,P2,Frob,flagTC,flagETCeffective,0,flagETCnonVanish,listeRpzGx,W2plusP0);
		J=resETC[1];
		listeRpzLambdaGAdic=[listeRpzLambdaGAdic,resETC[2]];
		flaglisteLambdaGAdic=[flaglisteLambdaGAdic,resETC[3]];
		);
	print("flag listeLambdaGAdic = ",flaglisteLambdaGAdic);
	return([J,listeRpzLambdaGAdic,flaglisteLambdaGAdic]);		
}





\\************************************** THETA FRACTION & HYPERPLANE **********************************************

\\flagTFEonephase=1; \\flag=0 IF we don't want to try TFE in one phase

InitVTFE(WK0,dataRRspace,flagETCeffective)=
{
	my(VTFE,VK0t);
	if(flagETCeffective==0,VTFE=WK0,
		VK0t=dataRRspace[5];
		VTFE=[WK0,VK0t];
		);
	return(VTFE);
}

InitHTFE(Omega,dataRRspace,flagETCeffective)=
{
	my(HK0t,HTFE);
	if(flagETCeffective==0,HTFE=Omega,
		HK0t=dataRRspace[4];
		HTFE=[Omega,HK0t];
		);
	return(HTFE);
}

InitThetaFractionOnePhase(J,listeFieldArchiOTC,listeRpzTFECoefBasis,listeOrbiteArchiAlpha,correspondance,HTFE,RR3D0,listeRpzGx,listeNum,flagTFEonephase,flagTC,flagETCeffective)=
{
	my(lCoef,nbO,listeTFE,failorbit,listeThetaFraction,nbE,failorbitEven,listeThetaFractionEven,listeTFEEven);
	if(flagETCeffective==0,lCoef=listeRpzTFECoefBasis,
                lCoef=listeRpzTFECoefBasis[1]);
	if(flagTFEonephase==0,
		if(flagTC,nbO=#listeFieldArchiOTC[1],nbO=#listeFieldArchiOTC);
        	failorbit=vector(nbO,i,i);
        	listeThetaFraction=vector(nbO,i,0);
		);

	if(flagTFEonephase,
		if(flagTC,[listeTFE,failorbit]=InitlisteTFE(J,listeFieldArchiOTC[1],lCoef,listeOrbiteArchiAlpha[1],correspondance[1],HTFE[1]),
			[listeTFE,failorbit]=InitlisteTFE(J,listeFieldArchiOTC,lCoef,listeOrbiteArchiAlpha,correspondance,HTFE));
		if(flagTC,listeThetaFraction = listeThetaPolynomialCore(J,listeRpzGx[1],listeNum[1],listeOrbiteArchiAlpha[1],listeFieldArchiOTC[1],correspondance[1],listeTFE,failorbit,RR3D0),
        		listeThetaFraction = listeThetaPolynomialCore(J,listeRpzGx,listeNum,listeOrbiteArchiAlpha,listeFieldArchiOTC,correspondance,listeTFE,failorbit,RR3D0));
		);

	if(flagTC&&flagETCeffective,
		[listeTFEEven,failorbitEven]=InitlisteTFE(J,listeFieldArchiOTC[2],listeRpzTFECoefBasis[2],listeOrbiteArchiAlpha[2],correspondance[2],HTFE[2]);
		listeThetaFractionEven = listeThetaPolynomialCore(J,listeRpzGx[2],listeNum[2],listeOrbiteArchiAlpha[2],listeFieldArchiOTC[2],correspondance[2],listeTFEEven,failorbitEven,RR3D0);
		);
	
	if(flagTC,
		nbE=#listeFieldArchiOTC[2];
		if(flagETCeffective==0,failorbitEven=vector(nbE,i,i));
		if(flagETCeffective==0,listeThetaFractionEven=vector(nbE,i,0));
		failorbit=[failorbit,failorbitEven];
		listeThetaFraction=[listeThetaFraction,listeThetaFractionEven];
		);
	return([listeThetaFraction,failorbit]);
}

/*
MakeTFEffectiveMain(listeThetaFractionInit,listeLambdaArchi,PolCurve,failorbit)=
{
                return(InitlisteTFEffective(PolCurve,listeThetaFractionInit,listeLambdaArchi,failorbit)[1]);
}
*/
/*
InitThetaFractionCore(J,listeRpzGx,listeNum,listeOrbiteArchiAlpha,listeFieldArchiOTC,correspondance,listeTFEInit,RR3D0,failorbit,flagTC)=
{
	my(listeThetaFraction);
	if(flagTC,listeThetaFraction = listeThetaPolynomialCore(J,listeRpzGx[1],listeNum[1],listeOrbiteArchiAlpha[1],listeFieldArchiOTC[1],correspondance[1],listeTFEInit[1],failorbit[1],RR3D0),
		listeThetaFraction=listeThetaPolynomialCore(J,listeRpzGx,listeNum,listeOrbiteArchiAlpha,listeFieldArchiOTC,correspondance,listeTFEInit,failorbit,RR3D0);
                );
	return(listeThetaFraction);
}
*/

InitThetaFractionMain(J,listeFieldArchiOTC,correspondance,listeOrbiteArchiAlpha,listeRpzGx,listeNum,listeTFEInit,RR3D0,failorbit,flagTFEonephase,flagTC)=
\\If TFE in one phase failed or flagTFEonephase=0 or flagTC=1
{
	my(flagOldSchool,listeThetaFraction);
        
	\\listeOrbiteLambdaGAdic=InitlisteOrbiteLambdaGAdic(J,listeFieldArchiOTC,correspondance,listeRpzLambdaGAdic);
	\\listeLambdaArchi=InitlisteLambdaArchi(J,listeOrbiteLambdaGAdic,correspondance,listeRpzOrbiteOTC,listeOrbiteArchiAlpha,listeFieldArchiOTC,WArfGAdic,failorbit);
	listeThetaFraction=listeTFEInit;

	\\******************* OTC *******************
	if(flagTC==0,flagOldSchool=(failorbit!=[])||(flagTFEonephase==0),
		flagOldSchool=(failorbit[1]!=[])||(flagTFEonephase==0));
	
	if(flagOldSchool,
		\\depend on flagTC
		\\listeThetaFraction=InitThetaFractionCore(J,listeRpzGx,listeNum,listeOrbiteArchiAlpha,listeFieldArchiOTC,correspondance,listeTFEInit,RR3D0,failorbit,flagTC);
		print("\n","COMPUTATION Odd Theta Fraction in the OLD SCHOOL way");
		if(flagTC,listeThetaFraction[1] = listeThetaPolynomialCore(J,listeRpzGx[1],listeNum[1],listeOrbiteArchiAlpha[1],listeFieldArchiOTC[1],correspondance[1],listeTFEInit[1],failorbit[1],RR3D0),
                	listeThetaFraction=listeThetaPolynomialCore(J,listeRpzGx,listeNum,listeOrbiteArchiAlpha,listeFieldArchiOTC,correspondance,listeTFEInit,failorbit,RR3D0);
                	);
		);

	\\****************** ETC *******************
	if(flagTC,
		flagOldSchool=(failorbit[2]!=[])||(flagTFEonephase==0);
		if(flagOldSchool,
			print("\n","COMPUTATION Even Theta Fraction in the OLD SCHOOL way");
			listeThetaFraction[2] = listeThetaPolynomialCore(J,listeRpzGx[2],listeNum[2],listeOrbiteArchiAlpha[2],listeFieldArchiOTC[2],correspondance[2],listeTFEInit[2],failorbit[2],RR3D0);
			\\listeThetaFractionEven=InitThetaFractionCore(J,listeRpzGx[2],listeNum[2],listeOrbiteArchiAlpha[2],listeFieldArchiOTC[2],correspondance[2],listeTFEInit[2],RR3D0,failorbit[2],0);
			);
		);

	return(listeThetaFraction);
}

InitFlagTFMain(listeThetaFraction,flagTC)=
{
	my(flagTF);
	if(flagTC==0,flagTF=InitFlagTF(listeThetaFraction),
		flagTF=InitFlagTF(listeThetaFraction[1])&&InitFlagTF(listeThetaFraction[2]));
	return(flagTF);
}



InitlisteTHMain(listeThetaFraction,HTFE,listeFieldArchiOTC,listeIndexTFE,flagTC)=
{
        my(listeTH);
	print("\n","Theta Hyperplanes");
        if(flagTC==0,listeTH=InitlisteTH(listeThetaFraction,HTFE,listeFieldArchiOTC,listeIndexTFE),
		listeTH=InitlisteTH(listeThetaFraction[1],HTFE[1],listeFieldArchiOTC[1],listeIndexTFE[1]));
	print("end");
        return(listeTH);
}


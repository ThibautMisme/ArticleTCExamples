\\****************************************************************** 
\\*                     MAIN: Computing TC                         *
\\******************************************************************


\\********************************************* INITIALISATION ***************************************************
	gettime();
	getwalltime();
        if((flagmethod=="TF")||(flagmethod=="TFwithCL"),flagEffective=1);
	if((flagmethodTC=="TF")||(flagmethodTC=="TFwithCL"),flagETCeffective=1);
	if(flagTH,flagEffective=1);
	if(flagEffective==0,flagTFEonephase=0;flagTH=0;);
	if(flagpopti==0,flagemptytors=0);
	if(flagTC==0,flagETCeffective=0;flagETCnonVanish=0;);
	if(flagETCeffective==0,flagETCnonVanish=0);

	print(" name curve = ",name);
        if(default(debug)>1,PrintInput(flagTC,flagmethod,flagTFEonephase,flagTH,flagEffective,flagmethodChiJ2,flagmethodTC));
	print("emax = ",emax);
	if(#variables(PolCurveInit)==3,PolCurveInit=subst(PolCurveInit,z,1));
        if(flagMC,[PolCurveInit,listPoint,P0Init,P1Init,P2Init]=MCIntervertxy(PolCurveInit,listPoint,P0Init,P1Init,P2Init));
	[PolCurve,P0,change,P1,P2]=MakeP0Infinite(PolCurveInit,P0Init,flagP0,P1Init,P2Init);
        if(Badprime==[],Badprime=CurveBadPrime(PolCurve,g));
        MakdisiConstFF=InitMakdisiConstFF(PolCurve,pInit,Badprime,flagemptytors,flagpopti,boundp);
        [p,aMak,qTor,CardJp,CardJq,LpFrob,emptytors]=MakdisiConstFF;
        [Omega,w,CCurve]=InitCurveConst(PolCurve);
J=InitMakdisiGAdic(PolCurve,Omega,w,CCurve,MakdisiConstFF,emax,P0,kind,Auts,flagPointJFp);
        d0=getd0(J);
        dw=getdw(J);
        g=getg(J);

[J,WK0tGAdic,WK0]=InitWK0tGAdic(J,Omega,CCurve,P0);
        WK0tFF=RedPointJ(J,WK0tGAdic,1);

\\******************************************** COMPUTING OVER Fq ****************************************************
        JFF=InitMakdisiFF(J);
        CCurvep=crvinit(Mod(PolCurve,getpMak(J)),'t,'a);
SBasisJ2FF=InitSympBasisJ2FF(JFF);
	if(flagtestbasisSymp,B=SBasisJ2FF;print("test Symplectic Basis Nicolas: ",lift(InitWeilMatrixNicolas(JFF,B)),"\B"););\\only test: not neccessary
WB0tFF=InitWB0tFF(JFF,CardJp,CardJq,aMak,WK0tFF,SBasisJ2FF,emptytors);

\\******************************************** LIFT TOWARD Qq ***********************************************************
[listeOrbiteOdd,listOTCNt,Arfdata,B0data,Frob,flagArf]=DataOTCLocal(JFF,CCurvep,P0,SBasisJ2FF,WB0tFF,flagTC);
[J,WB0tGAdic,SBasisJ2GAdic,PicTorsFrobGendata]=DataOTCGAdic(J,JFF,SBasisJ2FF,Frob,WB0tFF,WK0tGAdic,CardJq,CardJp);

listeRpzOrbiteOTC=InitRpzOrbiteOTC(J,listeOrbiteOdd,SBasisJ2GAdic,WB0tGAdic,flagTC);

\\*************************************************** INTERPOLATION ********************************************************
\\ChiOTC by alpha method
if((flagmethod=="Classic")||(flagmethod=="Nicolas"),[J,ChiOTC,listeOrbiteAlpha,flagChiOTC]=InitChiOTCAlphaMain(J,listeRpzOrbiteOTC,listeOrbiteOdd,flagmethod,flagTC,CCurve,P0,P1,P2,Frob));

\\************************************************ THETA FRACTIONS ***********************************************************
	print("\n","Theta Fractions: Denominator");
	\\IF problem: J=PrepChi(J,CCurve,P0,P1,P2);
	if(flagAlphaGlobal==0,W2plusP0=InitW2plusP0(J,CCurve,P0),W2plusP0=0);
listeRpzGx=InitlisteRpzGx(J,listeRpzOrbiteOTC,flagTC,W2plusP0);

	dataRRspace=InitdataRRspace(J,Omega,CCurve,P0);
	print("\n","Theta Fractions: Numerator");
listeNum=InitlisteNum(J,listeRpzOrbiteOTC,dataRRspace,flagTC);
	\\dataRRspace=[HK0tD0,K_K0tD0,VK0tD0,HK0t,VK0t,DK0t,RR3D0];	RR3D0=[H3D0,W3];
	RR3D0=dataRRspace[7];

	print("\n","Theta Fractions: WArf");
	WArfGAdic=InitWArfGAdic(J,CCurve,P0,flagETCeffective);
	flagAlphaGlobal=CheckflagAlpha(listeRpzGx,flagTC);
	if((flagAlphaGlobal==0),WArfGAdic=WArfGAdicflagAlpha(J,P0,WArfGAdic,flagETCeffective,flagETCnonVanish));
	dimArf=[0,0,0,0];

\\Lambda if flagEffective=1;
	print("\n","Theta Fractions: Lambda"); 
if(flagEffective,[J,listeRpzLambdaGAdic,flaglisteLambdaGAdic]=InitLambdaGAdicMain(J,listeRpzOrbiteOTC,WArfGAdic,dimArf,CCurve,P0,P1,P2,Frob,flagTC,flagETCeffective,flagETCnonVanish,listeRpzGx,W2plusP0));        

\\Lambda Coef if flagTFEonephase=1 (and flagEffective=1)
	if(flagTFEonephase&&flagEffective,VTFE=InitVTFE(WK0,dataRRspace,flagETCeffective));
if(flagTFEonephase&&flagEffective,listeRpzTFECoefBasis=InitlisteRpzTFECoefBasis(J,listeRpzLambdaGAdic,listeRpzGx,listeNum,VTFE,flagTC,flagETCeffective));

	print("\n","Etale Polynomial Chi");
\\Chi byTFE (only if flagEffective=1 AND flagTFEonephase=1)
if((flagmethod=="TF")||(flagmethod=="TFwithCL"),[ChiOTC,listeOrbiteAlpha,flagChiOTC]=InitChiOTCbyTFEMain(J,listeRpzTFECoefBasis,listeOrbiteOdd,flagmethod,flagTC,flagETCeffective,1));

\\if Chi failed a first time
if(flagChiOTC==0,print("we try ChiOTC by second method: ",flagmethodbis);[J,ChiOTC,listeOrbiteAlpha,flagChiOTC]=InitChiOTCbyTFEMix(J,listeRpzOrbiteOTC,listeRpzTFECoefBasis,listeOrbiteOdd,flagmethodbis,flagTC,CCurve,P0,P1,P2,Frob,flagETCeffective));

	flagChiETC=0;
\\ChiETC by TFE
if(flagTC&&((flagmethodTC=="TF")||(flagmethodTC=="TFwithCL")),[ChiETC,listeOrbiteAlphaEven,flagChiETC]=InitChiETCbyTFEMain(J,listeRpzTFECoefBasis,listeOrbiteOdd,flagmethodTC,flagETCeffective));
\\ChiETC by Alpha if flagTC=1 or failed the first time
	if(flagTC&&(flagChiETC==0),print("we try ChiETC by second method: ",flagmethodTCbis);flagmethodTC=flagmethodTCbis;);
if(flagTC&&(flagChiETC==0),[J,ChiETC,listeOrbiteAlphaEven,flagChiETC]=InitChiETCAlpha(J,listeRpzOrbiteOTC,listeOrbiteOdd,flagmethodTC,CCurve,P0,P1,P2,Frob));
	if(flagTC,[ChiOTC,listeOrbiteAlpha,flagChiOTC]=UpdateChiOTC(ChiOTC,ChiETC,listeOrbiteAlpha,listeOrbiteAlphaEven,flagChiOTC,flagChiETC));


	print("\n","Theta Fractions");
	if(flagTC&&(flagChiETC==0),flagTC=0;flagETCeffective=0;listeNum=listeNum[1];listeRpzGx=listeRpzGx[1];);

[listeFieldArchiOTC,listeOrbiteArchiAlpha,correspondance,listsmallindex,flagThreshold] = InitArchidata(J,ChiOTC,listeOrbiteAlpha,flagTC,thresholdOrbit);

	if(flagThreshold,[listeFieldArchiOTC,listeOrbiteArchiAlpha,correspondance]=ThresholdOrbitClean(listsmallindex,listeFieldArchiOTC,listeOrbiteArchiAlpha,correspondance,thresholdOrbit,flagTC));

\\TF if flagTFEonephase=1. If it is a success, we are done (with the OTC) + the TF (Odd) are Effective. We set up variables otherwise
	HTFE=InitHTFE(Omega,dataRRspace,flagETCeffective);
[listeThetaFraction,failorbit]=InitThetaFractionOnePhase(J,listeFieldArchiOTC,listeRpzTFECoefBasis,listeOrbiteArchiAlpha,correspondance,HTFE,RR3D0,listeRpzGx,listeNum,flagTFEonephase,flagTC,flagETCeffective);

\\TF if flagTFEonephase=0 or TFE one phase failed (or flagTC=1)
listeThetaFraction=InitThetaFractionMain(J,listeFieldArchiOTC,correspondance,listeOrbiteArchiAlpha,listeRpzGx,listeNum,listeThetaFraction,RR3D0,failorbit,flagTFEonephase,flagTC);

\\We make the Theta Fractions Effective when we can
if(flagEffective,[listeThetaFraction,listeIndexTFE]=MakeTFEffectiveMain(J,PolCurve,listeThetaFraction,listeFieldArchiOTC,correspondance,listeRpzLambdaGAdic,listeOrbiteArchiAlpha,failorbit,flagTC,flagETCeffective));
	
	flagThetaFraction=InitFlagTFMain(listeThetaFraction,flagTC);
	print("Theta Fractions: end");

\\**************************************************** THETA HYPERPLANES ***************************************************************
if(flagEffective&&flagTH,listeTH=InitlisteTHMain(listeThetaFraction,HTFE,listeFieldArchiOTC,listeIndexTFE,flagTC));


\\******************************************* END ********************************************
        print("\n","total time (CPU): ",strtime(getabstime()));
	timecomputation=gettime();

	if(flagTC==0,chiOTC=ChiOTC[1],chiOTC=ChiOTC);
	if(flagTC,chiETC=ChiOTC[2]);
        g=getg(J);
	listTF=InitlistTF(listeThetaFraction,flagTC,flagThetaFraction);
	listFieldTC=listeFieldArchiOTC;

	print("\n","Goodbye ;)");

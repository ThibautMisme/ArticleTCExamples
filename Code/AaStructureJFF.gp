\\ \r AaTechnique.gp
\\ \r AaMakdisiMachinery.gp
\\ \r AaPairing.gp

\\print("Optimiser le calcul des t et ti ?");
\\print("faire calcul de la structure de J(Fp) par une représentation de J(Fp) au lieu de celle de J(Fq) ? c'est algèbre linéaire de matrices de meme taille mais coef sont des nombres de Fp au lieu d'etre des polynomes de Fp"); Comment choper la représentation puisqu'on a pas de structure d'anneau et qu'on peut pas utiliser les matrices de multplications? On perd aussi avantage du Pairing de Weil.

RandomPointEvenJFF(JFF,listeRandom)=
\\[uN,char]=listeRandom
\\uN=u*N where N=|J(Fq)odd|or|J(Fp)odd| and u*N +v*|J(Fp/q)even|=1
\\char="Jq" or "Jp" 
{
	my(P,uN,char);
	print("computation new point");
	[uN,char]=listeRandom;
	if(char=="Jq",
		P=RandomPointJ(JFF););
	if(char=="Jp",
		P=RandomPointJFp(JFF););
	P=MulInt(JFF,P,uN);\\projection on the even part
	return(P);
}

OrderPointEvenJFF(JFF,h,logm)=
{
	my(P,listeth,flag,logorder,order);
	P=h;
	listeth=vector(logm,i,0);\\latter for finding relations
	flag=(CheckZero(JFF,P)==0);
	logorder=0;
	while(flag,
		logorder=logorder+1;
		listeth[logorder]=P;
		if(logorder>logm,
			print("WARNING: PB INITITIATION B0: POINT D'ORDRE TROP GRAND");
			); 
		P=SumJ(JFF,P,P);\\P=2*P
                if(CheckZero(JFF,P),
			flag=0;
			);
		);
	\\order=2^logorder;
	return([logorder,listeth]);
}


CheckRelations(JFF,h,xh,G,listen)=
{
	my(P,ni,bool);
	P=MulInt(JFF,h,2^xh);
	for(i=1,#listen,
		ni=listen[i];
		if(ni>0,
			P=SumJ(JFF,P,MulInt(JFF,G[i][1],ni));
			);		
		);
	bool=CheckZero(JFF,P);
	return(bool);
}

Divlisten(listen)=
{
	my(NEWlisten,flagZERO,flagDIV,ni,flagy);
	NEWlisten=vector(#listen,i,0);
	flagZERO=1;\\all the n_i are 0 ?
	flagDIV=0;\\there exists i s.t. n_i is not divisible by 2 ?
	for(i=1,#listen,
		ni=listen[i];
		if(ni!=0,
			flagZERO=0;
			if((ni%2)!=0,
				flagDIV=1;
				);
			if((ni%2)==0,
				NEWlisten[i]=ni/2;
				);
			);
		);
	flagy=1;
	if( (flagZERO==0)&&(flagDIV==1),
		flagy=0;
		NEWlisten=listen;\\in that case yh=xh and we want to keep listen untouched
		);
	return([flagy,NEWlisten]);
}
				

Recurset(JFF,G,t,listen)=
{
	my(NEWt,ni);
	print("WAENING: passer pointeur");
	NEWt=t;
	for(i=1,#listen,
		ni=listen[i];
		if(ni!=0,
			if(NEWt!=0,
				NEWt=SumJ(JFF,NEWt,MulInt(JFF,G[i][1],ni)););
			if(NEWt==0,
				print("error Recurset: t=0")
				);
			););
	return(NEWt);
}

Recurselisten(listen,G,listeLambda)=
{
	return(vector(#listen,i, ( listen[i]+2^(G[i][2]-1)*listeLambda[i] )%(2^G[i][2]) ));
}

FindLambda(JFF,t,listetG,SBasisJ2FF)=
\\FindLambda(JFF,t,listetG,SBasisJ2,Weildata,CCurve);
{
	my(dim,tvec,tGMat,tGj,Lambda,flag);
	
	print("t torsion ? ",CheckTorsion(JFF,t)," t nul ? ",CheckZero(JFF,t));

	dim=2*getg(JFF);
	\\tvec=vector(dim,i,WeilPairingFF(JFF,CCurve,t,SBasisJ2FF[i],Weildata))~;
	tvec=vector(dim,i,WeilPairingFFNicolas(JFF,t,SBasisJ2FF[i]))~;
	tGMat=matrix(dim,#listetG,i,j,Mod(0,2));
	for(j=1,#listetG,
		tGj=listetG[j];
		if(tGj!=0,
			for(i=1,dim,
				\\WeilPairingFF
				tGMat[i,j]=WeilPairingFFNicolas(JFF,tGj,SBasisJ2FF[i]);
				);
			);
		);
	Lambda=matinverseimage(tGMat,tvec);
	flag=(length(Lambda)!=0);print("Lambda found ? ",flag);
	Lambda=lift(Lambda)~;
	return([flag,Lambda]);
}

FindNewhLoop(JFF,logm,G,BasisJ2,hinit,xhinit,listethinit,listetG)=
\\FindNewhLoop(JFF,logm,G,BasisJ2,hinit,xhinit,listethinit,listetG,Weildata,CCurve)=
{
	my(h,xh,listeth,t,listen,listenDIV,flag,yh,listeLambda);
	\\Initialisation new h:
	h=hinit;
	xh=xhinit;
	listeth=listethinit;
	if(xh==1,
		t=h;);
	if(xh!=1,
		t=listeth[xh];);
        listen=vector(#G,i,0);
        flag=1;

        \\Find y:
        while(flag,
	
		print("\n");
		print("xh = ",xh);
		print("G = ",vector(#G,i,G[i][2]));
		print("listen = ",listen);
		print("Check Relations : ",CheckRelations(JFF,h,xh,G,listen));		
		

		listenDIV=listen;\\storage if division by 2 and Lambda not found
        	[flag,listen]=Divlisten(listen);\\check the divisibility of the n_i by 2
                print("listen = ",listen);
		if(flag==0,
                	yh=xh;);
                if(flag,
			t=Recurset(JFF,G,t,listen);
                        \\[flag,listeLambda]=FindLambda(JFF,t,listetG,BasisJ2,Weildata,CCurve);
                        [flag,listeLambda]=FindLambda(JFF,t,listetG,BasisJ2);
			if(flag==0,\\Lambda_i not found
                        	yh=xh;
				listen=listenDIV;
				print("listen = ",listen);
                                );
                        if(flag,
                        	xh=xh-1;
                                if(xh==0,
                                	yh=xh;
                                        flag=0;
					listen=Recurselisten(listen,G,listeLambda);\\for DLJacEven
					print("listen = ",listen);
                                        );
                                 if(flag,
                                 	listen=Recurselisten(listen,G,listeLambda);
					print("listen = ",listen);
                                        if(xh==1,
						t=h;);
					if(xh!=1,
						t=listeth[xh];);
                                        );
                        	);
                   	);
		);
	if(yh!=0,
		print("WARNING: CheckRelation non nécessaire");
        	if(CheckRelations(JFF,h,yh,G,listen)==0,
        		print("BUG: Calcul Structure Jeven: mauvaise relation"););
		);

	return([h,yh,listen]);
}

InitlistetG(JFF,G)=
{
	my(listetG,gam);
	listetG=vector(#G,i,0);
	for(i=1,#G,
		gam=2^(G[i][2]-1);
		if(gam!=0,
			listetG[i]=MulInt(JFF,G[i][1],gam);
			);
		if(gam==0,
			listetG[i]=0;
			);
		);
        return(listetG);
}	

FindNewh(JFF,logm,G,BasisJ2,listeRandom)=
{
	my(listetG,yh,h,listen);

	listetG=InitlistetG(JFF,G);

	yh=0;
	while(yh==0,
		\\Initialisation new h:
        	xh=0;
        	while(xh==0,\\out of the FindNewhLoop so we can use it for DLJacEven
                	h=RandomPointEvenJFF(JFF,listeRandom);
                	[xh,listeth]=OrderPointEvenJFF(JFF,h,logm);
                	);
		[h,yh,listen]=FindNewhLoop(JFF,logm,G,BasisJ2,h,xh,listeth,listetG);
		);
	return([h,yh,listen]);
}

normaliseG(JFF,Ginitial,listen)=
{
	my(G,listev,ni,vi,ui);
	G=Ginitial;
	listev=vector(#listen,i,0);\\WARNING: now listev represents logarithm and not coefficients of a relation as listen used to do it.
	print("WARNING: A AMELIORER, PASSER DIRECTEMENT LE POINTEUR. cf évite copie. Changer gi en ui*gi est couteux. AMELIORATION???");
	for(i=1,#listen,
		ni=listen[i];
		if(ni==0,
			listev[i]=-1;);
		if(ni!=0,
			vi=DLint(ni,2);
			ui=ni/(2^vi);
			G[i]=[MulInt(JFF,G[i][1],ui),G[i][2]];
			listev[i]=vi;
			);
		);
	return([G,listev]);
}

CheckhGood(JFF,h,yh)=
{
	my(m);
	m=2^(yh-1);
	return(CheckZero(JFF,MulInt(JFF,h,2^yh))&&(CheckZero(JFF,MulInt(JFF,h,m))==0));
}
		
DiagonaliseCase1(JFF,Ginit,hinit,yh,listev)=
{
	my(G,h,vi,power);
	print("Diagonalise Case 1 :PASSER POINTEUR");
	G=Ginit;
	h=hinit;
	for(i=1,#listev,
		vi=listev[i];
		if(vi>-1,
			power=vi-yh;
			if(power==0,
				h=SumJ(JFF,h,G[i][1]););
			if(power!=0,
				h=SumJ(JFF,h,MulInt(JFF,G[i][1],2^power)););
			);
		);

	print("nouveau élément de G bon ?: ",CheckhGood(JFF,h,yh));

	G=concat(G,[[h,yh]]);
	return(G);
}

DiagonaliseCase2(JFF,Ginit,h,yh,listev)=
{
	my(G,mu,i0,gprime,vi,W0,hp,gamma0,yhp,listevp);
	G=Ginit;
	print("Diagonalise Case 2 :WARNING: pointeur");
	[mu,i0]=minv(listev);
	gprime=Ginit[i0][1];
	if(yh>mu,
		gprime=SumJ(JFF,gprime,MulInt(JFF,h,2^(yh-mu))););
	for(i=1,#G,
		if(i!=i0,
			vi=listev[i];
			if((vi>=0)&&(vi>mu),
				gprime=SumJ(JFF,gprime,MulInt(JFF,G[i][1],2^(vi-mu))););
			););
	W0=getW0(JFF);
	for(i=1,#G,
		if(i==i0,
			G[i]=[gprime,mu];);
		if(i!=i0,
			if(listev[i]>=0,
				G[i]=[AddFlipTM(JFF,Ginit[i][1],W0),Ginit[i][2]];);\\gi=-gi
			);
		);
	hp=AddFlipTM(JFF,h,W0);\\h=-h
	gamma0=Ginit[i0][2];
	yhp=yh+gamma0-mu;
	listevp=vector(#G,i,-1);
	for(i=1,#listevp,
		vi=listev[i];
		if((i!=i0)&&(vi>=0)&&((vi+gamma0-mu)<Ginit[i][2]),
			listevp[i]=vi+gamma0-mu;
			);
		);
	listevp[i0]=-1;

	print("G = ",vector(#G,i,G[i][2]));
	print("listevp = ",listevp);
	if(mu==0,\\then, g_i0 = 0
		G=concat(G[1..(i0-1)],G[(i0+1)..#G]);
		print("G = ",vector(#G,i,G[i][2]));
		listevp=concat(listevp[1..(i0-1)],listevp[(i0+1)..#listevp]);
		print("listevp = ",listevp);
		);

	return([G,hp,yhp,listevp]);
}
	
	

Diagonalise(JFF,Ginit,hinit,yhinit,listen)=
{
	my(G,listev,Case,vi,h,yh);
	\\make the relations better
	h=hinit;
	yh=yhinit;
	[G,listev]=normaliseG(JFF,Ginit,listen);
	print("listev = ",listev);
	
	Case=0;
	while(Case==0,
		Case=1;
		for(i=1,#listev,
			vi=listev[i];
			if((vi>=0)&&(yh>vi),
				Case=0;
				);
			);
	
		if(Case,\\Case 1
			G=DiagonaliseCase1(JFF,G,h,yh,listev);
			print("listev = ",listev);	
			);
		if(Case==0,\\Case 2
			[G,h,yh,listev]=DiagonaliseCase2(JFF,G,h,yh,listev);
			print("listev = ",listev);
			);
		);
	return(G);
}

InitG(JFF,logm,listeRandom)=
{
	my(gamma,g,listeh,G);
	print("WARNING: STOCKER LES PUISSANCES 2-ieme des G[i] ???");
        \\Initialisation G:
        gamma=0;
        while(gamma==0,
                g=RandomPointEvenJFF(JFF,listeRandom);
                [gamma,listeth]=OrderPointEvenJFF(JFF,g,logm);
                );
        G=[[g,gamma]];\\G=list of the basis. G[k]=[g,gamma] such that 2^gamma*g=0 is a relation

        /*listetG=[tg];\\WARNING: listetG is not the gathering of the listeth of the element g of G*/
	return(G);
}	

StructureJFF(JFF,m,listeRandom,BasisJ2)=
\\m=|J(Fp/q)even|
{
	my(logm,G,logmG,h,yh,listen);
	logm=logint(m,2);
	G=InitG(JFF,logm,listeRandom);
	logmG=G[1][2];

	\\print("BasisJ2 est elle bien symplectique ???");
	
	while(logmG!=logm,
		\\Find New h + Relations

		print("debut recherche New h");

		[h,yh,listen]=FindNewh(JFF,logm,G,BasisJ2,listeRandom);
			
		\\diagonalise new relation
		print("\n","debut diagonalise");
		G=Diagonalise(JFF,G,h,yh,listen);
		
		logmG=0;
		for(k=1,#G,
			logmG=logmG+G[k][2];
			);
		print("StructureFF: ecart de dimensions: ",logm-logmG);
		print("G : ",vector(#G,i,G[i][2]));
		);
	return(G);
}




DLJacEven(JFF,G,listetG,logm,BasisJ2,W)=
\\DLJacEven(JFF,CCurve,G,listetG,logm,BasisJ2,W,Weildata)=
\\m=|Group representend by G|
\\logm=log(m,2)
{
	my(xh,listeth,h,yh,listen);

	[xh,listeth]=OrderPointEvenJFF(JFF,W,logm);
	[h,yh,listen]=FindNewhLoop(JFF,logm,G,BasisJ2,W,xh,listeth,listetG);

	if(yh!=0,
		print("ERROR: DLJacEven: yW not equal to 0"););
	listen=vector(#listen,i,(-listen[i])%(2^G[i][2]));
	return(listen);
}








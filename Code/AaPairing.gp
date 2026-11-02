\\ \r AaTechnique.gp
\\ \r AaMakdisiMachinery.gp

\\print("Weil Pairing: ne pas calculer Delta1 ? puisqu'on prend quotient"); On prend celui de Nicolas de toute facon
\\print("calcul listeBranchP: comment faire avec changement variable?"); Nicolas
\\print("Optimisation: D0 autre chose que d0P0 ?"); Non trop compliqué: on a besoin de petit degré pour pouvoir créer les différents RR spaces
\\print("Optimisation: changer rpz des points avec un minimum de calculs"); Nicolas
\\print("Optimisation: stocker bonnes rpz et bons fP ?"); Nicolas

/*
ProofFunctionTors(J,W)=
\\for W s.t. [2]W~W0 (2-torsion point), returns f s.t. [2]W=W0+div(f) as a column of values
{
        my(W0,Wp,WZ,Kp,bool,f);
        W0=getW0(J);
        Wp=SumJ(J,W,W);
        Kp=WtoK(J,Wp);
        WZ=DivSubTM(J,Wp,W0,1,Kp);
        bool=length(WZ);
        f=0;
        if(bool!=1,
                print("ERROR: W is not torsion"););
        if(bool==1,
                f=WZ[,1];);
        return(f);
}
*/
/*
ProofFunctionTorsFF(JFF,W)=
\\for W s.t. [2]W~W0 (2-torsion point), returns f s.t. [2]W=W0+div(f) as a column of values
{
	my(W0,Wp,WZ,Kp,WZp,bool,f);
        W0=getW0(JFF);
	Wp=SumJ(JFF,W,W);
        WZ=DivSubDimFF(JFF,Wp,W0);

	\\WZ=DivSubDimRandomFF(JFF,Wp);	

	Kp=WtoK(JFF,Wp);
	WZp=DivSubTM(JFF,Wp,W0,1,Kp);
	print("DivSub pareil: ",(WZ)==(WZp));
        bool=length(WZ);
	f=0;
	if(bool!=1,
		print("ERROR: W is not torsion"););
	if(bool==1,
		f=WZ[,1];);
	return(f);
}
*/
ProofFunctionTorsFF2(JFF,W)=
{
        my(d0,g,d,Wp,W2,WZ,bool,f);
        d0=getd0(JFF);
        g=getg(JFF);
        d=2*d0-g+1;
        Wp=DivAddTM(JFF,W,W,d);\\4D0-2D
        W2=getWn(JFF,2);
        WZ=DivSubDimFF(JFF,Wp,W2);\\2D0-2D
        bool=length(WZ);
        f=0;
        if(bool!=1,
                print("ERROR: W is not torsion"););
        if(bool==1,
                f=WZ[,1];);
        return(f);
}
/*
EvaluationSectionDivisorFFThibaut(JFF,W,f)=
\\return f(D) if W represent D and f is given as a column of values
{
	my(MV,Mats,d0,dwp,Matwp,Matt,A,K,Ks,Kt,d1,Delta1,d2,Delta2);
	MV=getWn(JFF,2);
	Mats=FindSupplementFF(JFF,W,MV);
	d0=getd0(JFF);
	dwp=4*getd0(JFF)+1-getg(JFF);
	Matwp=DivAddTM(JFF,W,getWn(JFF,3),dwp);
	Matt=FindSupplementV5FF(JFF,Matwp);
	A=concat([Mats,Matt,Matwp]);
	K=KtoWFF(JFF,A);
	Ks=K[1..d0,];
	Kt=K[(d0+1)..(2*d0),];
	d1=matdet(reduceFFMatrix(JFF,Ks));print("d1 = ",liftFF(d1));
	\\print("delta1 obtenu ",d1);
	if(d1==0,
		print("FAIL IN EvaluationSectionDivisor: d1"););
	if(d1=!0,
		Delta1=matdet(reduceFFMatrix(JFF,Kt))/d1;print("Delta1= ",liftFF(Delta1));
		if(Delta1==0,
			print("FAIL IN EvaluationSectionDivisor: Delta1"););
		);
	A=concat([MulPointWiseFF(JFF,f,Mats),Matt,Matwp]);
	K=KtoWFF(JFF,A);
        Ks=K[1..d0,];
        Kt=K[(d0+1)..(2*d0),];
        d2=matdet(reduceFFMatrix(JFF,Ks));print("d2 = ",liftFF(d2));
        if(d2==0,
                print("FAIL IN EvaluationSectionDivisor: d2",d2););
        if(d2=!0,
                Delta2=matdet(reduceFFMatrix(JFF,Kt))/d2;print("Delta2 = ",lift(Delta2));
                if(Delta2==0,
                        print("FAIL IN EvaluationSectionDivisor: Delta2",liftFF(Delta2)););
                );
	return(Delta2/Delta1);
}
*/
EvaluationSectionDivisorFF(JFF,W,f)={return(reduceFF(JFF,EvalDivisorNicolas(JFF,W,f)));}
/*
WeilPairingFF(JFF,P,Q)=
{
	my(fP,fQ,c,fPQ,fQP,Weil,res);
	fP=ProofFunctionTorsFF2(JFF,P);
	fQ=ProofFunctionTorsFF2(JFF,Q);
	\\c=RandomPolCol(JFF);
	\\fP=MulVectFF(JFF,fP,c);
	\\fP=MulVectFF(JFF,fQ,c);

	fPQ=EvaluationSectionDivisorFF(JFF,Q,fP); print("fPQ: ",liftFF(fPQ));
	fQP=EvaluationSectionDivisorFF(JFF,P,fQ); print("fQp: ",liftFF(fQP));
	Weil=fPQ/fQP;
	res=Weil;
	print("weil : ",liftFF(Weil));
	if(Weil==reduceFF(JFF,1),
		res=Mod(0,2););
	if(Weil==reduceFF(JFF,-1),
		res=Mod(1,2););
	return(res);
}

WeilRandom(JFF,P,Qinitial)=
{
        my(Q,fP,eva,fPQ,fQ,fQP);
        Q=Qinitial;
        fP=ProofFunctionTorsFF2(JFF,P);
        eva=liftFF(EvaluationSectionDivisorFF(JFF,Q,fP));
        while(eva==0,
                Q=SumJRandom(JFF,Q,EquivalentW0(JFF));
                eva=liftFF(EvaluationSectionDivisorFF(JFF,Q,fP));
                );
        fPQ=reduceFF(JFF,eva);
        fQ=ProofFunctionTorsFF2(JFF,Q);
        fQP=EvaluationSectionDivisorFF(JFF,P,fQ);
}
*/

domfP(J,CCurve,listP,f)=
\\flag=0 if dom=+oo
{
        my(varf,dom,SumVal,z,t,B,Val,fB,flag);

        \\print("obtenir listeBranch ? cf Changement de variable");
        varf=variables(JFF[4][2]);
        /*
        varf=variables(f);
        if(#varf==1,
                varf=[varf[1],y];);
        if(#varf==0,
                varf=[x,y]);
        */
        \\z=varlower("z");
        dom=1;
        SumVal=0;
	flag=1;
        \\for(k=1,#listeBranchP,
	for(k=1,#listP,
                \\B=listeBranchP[k];
		B=listP[k];
                \\print("Pt2Branch possible sur Fq ? gerer les variables");
                \\Val=BranchValuation(f,B[1],varf[1],varf[2]);\\print("Val = ",Val);
                Val=crvfnord(CCurve,f,B);
		Val=-Val;
                SumVal=SumVal+Val;
                \\fB=BranchEval(f,B[1],Val+1,varf[1],varf[2]);\\print("fB = ",fB);
                fB=crvfneval(CCurve,f,B);
		z=variable(fB);\\print("variable fB = ",z);
                if((fB!=+oo)&&flag,
			flag=1;
			fB=Pol((z^Val)*fB);
                	\\print("fB = ",liftFF(fB));
                	dom=dom*polcoef(fB,0);\\print("dom = ",liftFF(dom));
			);
		if(fB==+oo,
			flag=0;
			dom=+oo;
			);
                );
        if(SumVal!=2*getd0(J),
                print("error: proof function of degree too low");
                );
        return([flag,dom]);
}

/*
Inithdata(JFF,listeBranchP)=
{
        my(valh,c,hpol,ch,hWeil,WWeil);
        valh=0;
        while(valh!=getd0(JFF),
                c=vector(#getW0(JFF),i,random(30))~;
                hpol=JFF[4][1]*c;
                ch=PoltoMakdisiColFF(JFF,hpol);
                valh=BranchValuation(hpol,listeBranchP[1][1],x,y);
                valh=-valh;
                );
        hWeil=[hpol,ch];
        WWeil=mulPointWise(JFF,ch,getW0(JFF));
        return([hWeil,WWeil]);
}
*/

Inithdata(JFF,CCurve,P)=
\\P can be a point or a branch (if singular point)
{
        my(valh,c,hpol,ch,hWeil,WWeil);
        valh=0;
        while(valh!=getd0(JFF),
                c=vector(#getW0(JFF),i,random(30))~;
                hpol=JFF[4][1]*c;
                ch=PoltoMakdisiColFF(JFF,hpol);
                \\valh=BranchValuation(hpol,listeBranchP[1][1],x,y);
                valh=crvfnord(CCurve,hpol,P);
		valh=-valh;
                );
        hWeil=[hpol,ch];
        WWeil=mulPointWise(JFF,ch,getW0(JFF));
        return([hWeil,WWeil]);
}


WeilPairingFF(JFF,CCurve,Pinitial,Qinitial,Weildata)=
{
        my(Q,P,hWeil,WWeil,listP,h,hmod,ch,fP,fPg,fPQ,hP,WeilP,gP,fQ,fQpol,fQP,fQtD0,fQtD01,fQtD02,fQtD02pre,WeilQ,Weil);
        Q=Qinitial;
        P=Pinitial;
	\\[hWeil,WWeil,listeBranchP]=Weildata;
	[hWeil,WWeil,listP]=Weildata;
        [h,ch]=hWeil;

        fP=ProofFunctionTorsFF2(JFF,P);
        fPg=EvaluationSectionDivisorFF(JFF,WWeil,fP);
        fPQ=EvaluationSectionDivisorFF(JFF,Q,fP);
        hP=EvaluationSectionDivisorFF(JFF,P,ch);
        while((liftFF(fPg)==0)||(liftFF(fPQ)==0)||(liftFF(hP)==0),
                \\print("changement de rpz de P");
                P=SumJRandom(JFF,P,EquivalentW0(JFF));
                fP=ProofFunctionTorsFF2(JFF,P);
                fPg=EvaluationSectionDivisorFF(JFF,WWeil,fP);
                fPQ=EvaluationSectionDivisorFF(JFF,Q,fP);
                hP=EvaluationSectionDivisorFF(JFF,P,ch);
                );
        WeilP=fPQ/fPg;\\print("WeilP = ",liftFF(WeilP));

        gP=1/(hP^2);
        fQ=ProofFunctionTorsFF2(JFF,Q);
        fQpol=ColtoPolFFW2(JFF,fQ);
        fQP=EvaluationSectionDivisorFF(JFF,P,fQ);
        \\fQtD01=(domfP(JFF,listeBranchP,fQpol)/domfP(JFF,listeBranchP,h^2))^getd0(J);
	\\fQtD01=(domfP(JFF,CCurve,listP,fQpol)/domfP(JFF,CCurve,listP,h^2))^getd0(J);
	
	/*
	[flag1,fqp]=domfP(JFF,CCurve,listP,fQpol);
	[flag2,hp]=domfP(JFF,CCurve,listP,h^2);
	if(flag1&&flag2,fQtD01=(fqp/hp)^getd0(J),flag=0);
	*/

	\\if(BranchValuation(fQpol/h^2,listeBranchP[1][1],x,y)!=0,
	if(crvfnord(CCurve,fQpol/h^2,listP[1])!=0,
                print("error: fQpol/h^2 is not of valuation 0 over P0"););
        \\fQtD02pre=Pol(BranchEval(fQpol/h^2,listeBranchP[1][1],1,x,y));\\print("fQtD02pre = ",fQtD02pre);
        fQtD02pre=Pol(crvfneval(CCurve,fQpol/h^2,listP[1]));
	fQtD02=polcoef(fQtD02pre,0)^getd0(JFF);
        \\print("fQtD01==fQtD02 ", fQtD01==fQtD02);
        WeilQ=fQP*gP/fQtD02;\\print("WeilQ = ",liftFF(WeilQ));

        Weil=WeilP/WeilQ;\\print("Weil = ",liftFF(Weil));
	
	if(liftFF(Weil)==1,
                Weil=Mod(0,2);
                );
        if(liftFF(Weil)==(getpMak(JFF)-1),
                Weil=Mod(1,2);
                );
        return(Weil);
}

WeilPairingFFNicolas(JFF,P,Q)=
{
	my();
	weil=picweil(JFF,2,P,Q);
	if(weil==1,weil=Mod(0,2),
		weil=Mod(1,2);
		);
	return(weil);
}	

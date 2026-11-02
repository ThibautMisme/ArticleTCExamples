install("FpX_root_order_bound","GG")
\\ \r AaTechnique.gp

/* NOTATIONS
CCurve is a curve provided by the function CrvInit in Algcurves.gp

*/

InitLpFrob(CCurveFF)=
\\the curve provided has to p a curve in Fp
{
	return(crvcharpoly(CCurveFF));
}

InitaMak(LpFrob)=
{
	return(FpX_root_order_bound(LpFrob,2)[2]);
}

InitCardJac(LpFrob,aMak)=
{
	my(t,cardJp,cardJq);
	t=variable(LpFrob);
	cardJp=subst(LpFrob,t,1);
	cardJq=polresultant(LpFrob,t^aMak -1);
	return([cardJp,cardJq]);
}

InitMakdisiConstFF(PolCurve,pInit,Badprime,flagemptytors,flagpopti,boundp)=
{
	my(deg,Boundg,d0,nZ,flag,CC,g,p,pMak,index,listGoodPrime,CCurveFF,LpFrob,aMak,qTor,emptytors,CardJp,CardJq,MakdisiConstFF,goodindex,listemptytors,goodCardJp,goodCardJq);
	deg=poltotaldegree(Polcurve);
	Boundg=(deg-1)*(deg-2)/2;
	d0=2*Boundg+1;
	nZ=5*d0+1;
	flag=1;
	CC=crvinit(PolCurve,'t,'a);
	g=crvgenus(CC);
	\\p=max(pInit,11);\\Comment trouver la bonne borne afin d'assurer que C(Fp) soit suffisamment grand pour generer J(Fp) ? 
	p=pInit;

	if(flagpopti,print("we try to optimise p"));
	if(flagemptytors,print("we search for a prime with no rational 2-torsion"));

	\\CHECK good reduction mod p
	/*
	BadPrime=CurveBadPrime(PolCurve);print("Bad Primes Discriminant = ",BadPrime);
	index=0;
	for(k=1,#BadPrime,
		index=index+1;
		if(CheckSingular(PolCurve,g,BadPrime[index])==0,
			print("good p: ",BadPrime[index]);
			BadPrime=concat(BadPrime[1..(index-1)],BadPrime[(index+1)..#BadPrime]);
			index=index-1;
			);
		);
	print("Bad primes = ",BadPrime);
	\\[p,CCurveFF]=GoodRed(PolCurve,p,BadPrime);	
	*/

	pMak=0;
	listGoodPrime=[];
	while((#listGoodPrime<10)&&(p<boundp),
		CCurveFF=crvinit(PolCurve,'t,'a,Mod(1,p));
		flag=(CCurveFF==0);
		if(flag,print("CCurveFF failed, we increase p"),
			LpFrob=InitLpFrob(CCurveFF);\\charpoly of Frob_p acting on J[l]
			aMak=InitaMak(LpFrob);
			qTor=p^aMak;
		
			\\CHECK aMak
			flag=( (qTor/2)<(nZ+1) );
			
			if(flag,print("\n","WARNING: qTor too low to contain enough points to represent Makdisi setting. pMak increases:OLD pMak = ",p,"\n"),
				\\CHECK good reduction mod p
				flag=IsInList(Badprime,p);
				if(flag,print("WARNING: pMak is a prime of bad reduction. pMak increases:OLD pMak = ",p),
					pMak=p;\\this is a valid prime
					if(flagpopti,listGoodPrime=concat(listGoodPrime,[[pMak,CCurveFF,LpFrob,aMak,qTor]]),
						listGoodPrime=vector(10,i,0);
						);
					);
				);
			);
		p=nextprime(p+1);
		);

	if(flagpopti==0,
		emptytors=Mod(subst(LpFrob,variable(LpFrob),1),2);
		[CardJp,CardJq]=InitCardJac(LpFrob,aMak);
		);
	if(flagpopti,
		if(flagemptytors,
			listemptytors=[];
			for(k=1,#listGoodPrime,
				LpFrob=listGoodPrime[k][3];
				emptytors=Mod(subst(LpFrob,variable(LpFrob),1),2);\\emptytors==0 <-> J[l](Fp) NON zero \\emptytors==1 <-> J[l](Fp)=0
				if(emptytors==1,
					listemptytors=concat(listemptytors,[[k,emptytors]]));
				);
			if(#listemptytors!=0,print("we found p such that J[l](Fp)==0"),
                        	print("we couldn't find p such that J[l](Fp)==0");
                    		listemptytors=vector(#listGoodPrime,k,[k,Mod(0,2)]);
				);
			);
		if(flagemptytors==0,
			listemptytors=vector(#listGoodPrime,k,0);
			for(k=1,#listemptytors,
				LpFrob=listGoodPrime[k][3];
                                emptytors=Mod(subst(LpFrob,variable(LpFrob),1),2);
				listemptytors[k]=[k,emptytors];
                                );
			);
		goodindex=1;
		aMak=listGoodPrime[listemptytors[1][1]][4];
                LpFrob=listGoodPrime[listemptytors[1][1]][3];
                [CardJp,CardJq]=InitCardJac(LpFrob,aMak);
		goodCardJp=CardJp;
		goodCardJq=CardJq;
		for(k=1,#listemptytors,
			aMak=listGoodPrime[listemptytors[k][1]][4];
			LpFrob=listGoodPrime[listemptytors[k][1]][3];
			[CardJp,CardJq]=InitCardJac(LpFrob,aMak);
			print(" p = ",listGoodPrime[listemptytors[k][1]][1]," goodCardJq = ",goodCardJq," CardJq = ",CardJq);
			if(CardJq<goodCardJq,
				goodindex=k;
				goodCardJp=CardJp;
                		goodCardJq=CardJq;
				);
			);
		[pMak,CCurveFF,LpFrob,aMak,qTor]=listGoodPrime[listemptytors[goodindex][1]];
		CardJp=goodCardJp;
		CardJq=goodCardJq;
		emptytors=listemptytors[goodindex][2];
		);

	print("\n","final choices: pMak= ",pMak," qTor= ",qTor," aMak= ",aMak,"\n"); 
	print("cardinal of the Jacobian: over Fp = ",CardJp," over Fq = ",CardJq);
	if(emptytors!=0,print(" J[l](Fp)==0"),
                        print(" J[l](Fp)!=0");
                        );
	MakdisiConstFF=[pMak,aMak,qTor,CardJp,CardJq,LpFrob,emptytors];
	return(MakdisiConstFF);
}

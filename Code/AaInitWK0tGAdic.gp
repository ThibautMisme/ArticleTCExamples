\\ \r AaTechnique.gp
\\ \r AaMakdisiMachinery.gp

InitW2primepolGAdic(J,CCurve,P0)=
\\      W2prime rpz H°(D0 + (2g-2)P0)
/* return: a basis of H°(D0 + (2g-2)P0) over Qp in Makdisi setting */
\\IL FAUT QUE D0=d0P0
{
	my(d0,g,D2prime,W2primepol,W2prime);
	\\print("\n","WARNING: obtention de H°(d0P0 + (2g-2)P0) peut être optimisé","\n","IL FAUT QUE D0=d0P0","\n");
	d0=getd0(J);
	g=getg(J);
	D2prime=[P0,d0+2*g-2];
	W2primepol=RiemannRochTM(CCurve,D2prime);
	\\W2prime=PoltoMakdisi(J,W2primepol);
	return(W2primepol);
}

InitWK0GAdic(J,Omega)=
{
	my(WK0);
	WK0=PoltoMakdisi(J,Omega);\\print("rank WK0 = ",matrank(rpzJtoModMatrix(J,WK0))," = nb col = ",#WK0);
	return(WK0);
}

InitWK0tGAdic(Jinit,Omega,CCurve,P0)=
{
	my(J,WK0,dim,W2primepol,W2prime,K2prime,WK0t);
	J=Jinit;
	WK0=InitWK0GAdic(J,Omega);if(default(debug)>1,print("dim WK0 = ",#WK0));
	dim=getd0(J)-getg(J)+1;
	W2primepol=InitW2primepolGAdic(J,CCurve,P0);
	W2prime=PoltoMakdisi(J,W2primepol);\\print("rank W2prime = ",matrank(rpzJtoModMatrix(J,W2prime))," = nb col = ",#W2prime);
	K2prime=WtoKGAdic(J,W2prime);
	WK0t=DivSubTM(J,W2prime,WK0,dim,K2prime);
	
	\\if((flagP0==0)&&((#P0==2)||((#P0==3)&&(P0[3]==1))),
	if(CheckWinK(J,WK0t,getKn(J,2))==0,
		while((geteJ(J)>10)&&(CheckWinK(J,WK0t,getKn(J,2))==0),
			if(default(debug)>1,print("loss of precision (-5) because P0 is affine"));
			J[7]=geteJ(J)-5;
			);
		if(default(debug)>1,print("new precision: ",geteJ(J)));
		);
	if(default(debug)>1,print("WtoK good ??"));
	\\print("Check K2prime: ",CheckWinK(J,W2prime,K2prime));
	\\print("Check W2prime dans K2: ",CheckWinK(J,W2prime,getKn(J,2)));
	if(default(debug)>1,print("Check WK0t in V2: ",CheckWinK(J,WK0t,getKn(J,2))));
	\\It works because K0 is base point free: WK0t={v in W2prime | v*WK0 in W2prime}=H°(D0 + (2g-2)P0 - K0)
	if(default(debug)>1,print("WK0tGAdic : done \n"));
	return([J,WK0t,WK0]);	
}

/*
InitWK0tGAdicBourrin(J,Omega,CCurve,P0)=
{
	my(d,g,c,f,Df,index,DK0,flag,P0p,D,WK0pol,WK0);
	d=0;
	g=getg(J);
	while(d!=(2*g-2),
		d=0;
		c=vector(#Omega,i,1+random(30))~;
		f=Omega*c;
		\\Df=FnDivTM(CCurve,f);
		Df=crvfndiv(CCurve,f);
		index=1;
		while(Df[index,2]>=0,
			index=index+1;
			);
		DK0=Df[index..#Df[,1],];
		for(k=1,#DK0[,1],
			d=d+DK0[k,2];
			);
		d=-d;
		);
	print("d=",d);
	print("DK0 =",DK0);	
	flag=0;
	P0p=P0;
	if(#P0==2,
		P0p=[P0[1],P0[2],1];);
	for(k=1,#DK0[,1],
		if(CheckEgalitePointCoord(P0p,DK0[k,1]),
			flag=[DK0[k,2],k];
			);
		);
	if(flag!=0,
		d=d+flag[1];
		DK0=concat(DK0[1..(flag[2]-1),]~,DK0[(flag[2]+1)..(#DK0[,1]),]~)~;
		);

	D=concat([P0p,getd0(J)+d]~,DK0~)~;print("D=",D);
	WK0pol=RiemannRochTM(CCurve,D);print("WK0pol = ",WK0pol);
	WK0=PoltoMakdisi(J,WK0pol);
	return(WK0);
}	
*/



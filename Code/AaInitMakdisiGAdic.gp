\\install("PicInit","GGLLGGGGLGL")
\\install("SmoothRRdata","GGG");
\\ \r AaInitCurveConst


InitMakdisiGAdic(PolCurve,Omega,w,CCurve,MakdisiConstFF,emax,P0,kind,Auts,flagPointJFp)=
{
	my(pMak,aMak,LpFrob,emptytors,data,bad,J,d0,g);
	\\Version Nicolas
/*
	if(kind=="SmoothPlane",
		data=SmoothRRdata(PolCurveInit,pMak,P);
		);
	fdata=data[1];\\PolCurve after change of variables made by SmoothRRdata. Starting from there, no more changes of equation
*/	
	pMak=MakdisiConstFF[1];
	aMak=MakdisiConstFF[2];
	LpFrob=MakdisiConstFF[6];
	emptytors=MakdisiConstFF[7];
	\\Version D0=d0P0
	\\[Omega,PolCurve,w,CCurve]=InitCurveConst(PolCurve);
	
	\\bad=w;

	\\bad=w JUSTIFICATION: 
	\\We want PicInit to avoid using the poles of the elements of Omega as evaluation points since we will represent H°(K0)=Omega
	\\PicInit will pay attention to avoid the zeros of bad
	\\fi in Omega = wi/w with wi a POLYNOMIAL. So poles of fi < zeros of w UNION points at the infinity
	\\But PicInit don't take evaluation points at the infinity, so it is enough to tell him to avoid the zeros of w 	
	
	if(#P0==2,
		bad=bad*(x-P0[1]););
	if((#P0==3)&&(P0[3]==1),
		bad=bad*(x-P0[1]););
	\\we also force Z not to contain P0
	data=vector(5,i,0);
        data[1]=PolCurve;
        data[2]=Auts;\\Auts Automorphisme
	g=#Omega;
	data[3]=g;
	d0=2*g+1;
	data[4]=d0;\\print(" P0 = ",P0," d0 = ",d0);
	\\print("optimisation: tronquer base de H°(2D0) pour avoir base de H°(D0)");
	if(default(debug)>1,print("RR computation: V1, V2, VE1, VE2"));
	data[5]=[crvRR(CCurve,Mat([P0,d0])) , crvRR(CCurve,Mat([P0,2*d0])) , [InitWE(0,CCurve,d0,g,P0,P1,0),InitWE(0,CCurve,d0,g,P0,P2,0)] ];\\3ième élément est info pour E1,E2, cf fonction alpha

	if(default(debug)>1,print("RRspace: H0(D0), H0(2D0) ",data[5][1..2]));

	bad=denominator(data[5]);
	\\mfact=factor(bad);for(k=1,#mfact[,1],if(mfact[k,1]==x,bad=bad/(mfact[k,1]^mfact[k,2])));
	if(default(debug)>1,print("bad = ",bad));

	if(default(debug)>1,print("PicInit beginning"));
	if(default(debug)>1,print("if error: 'PicInit: impossible inverse in Flxq_inv: 0', look at the setting of bad in AaInitMakdisiGAdic.gp "));	
	\\picinit(f,Auts,g,d0,L,badlocus,p,a,e,Lp):
	if(flagPointJFp&&(emptytors==0),\\we add points of J(Fp) in Z
		if(default(debug)>1,print("PicInit Thibaut"));
		J=picinit(data[1],Auts,data[3],data[4],data[5],bad,pMak,aMak,emax,LpFrob,1);
		);
	if((flagPointJFp==0)||emptytors,
		if(default(debug)>1,print("Classic picinit"));
		if(default(debug)>1,print("flagPointJFp = ",flagPointJFp));
		J=picinit(PolCurve,Auts,g,d0,data[5],bad,pMak,aMak,emax,LpFrob,0);
		);

	if(default(debug)>1,print("PicInit end"));
	return(J);
}

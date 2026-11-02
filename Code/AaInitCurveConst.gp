\\ \r AaTechnique.gp
\\ \r AaMakdisiMachinery.gp

CheckP0(CCurve,P0)=
{
	return(crvfneval(CCurve,x,P0)==+oo);
}

MakeP0Infinite(PolCurveInit,P0Init,flagP0,P1Init,P2Init)=
\\make a change of variable (if flag==1) to force P0 to be on the chart z=0
\\change gives the coordinates of the map phi s.t. for f the homogeneous equation of the initial curve and f' the one of the new one,
\\we have:     f' ° phi = f
{
	my(CCurve,flag,PolCurve,P0,b,a,change,f,P1,P2);
	CCurve=crvinit(PolCurveInit,'t,'a);
	flag=flagP0;
	P1=P1Init;
	P2=P2Init;
	if(P1!=0,
		if(type(P1)=="t_INT",print("WARNING:P1 singular. P0 should be such that no change of variables is made"));
		if(#P1==2,
			P1=[P1[1],P1[2],1]));
	if(P2!=0,
		if(type(P2)=="t_INT",print("WARNING:P2 singular. P0 should be such that no change of variables is made"));
                if(#P2==2,
                        P2=[P2[1],P2[2],1]));
	if(#P0Init==3,P0=P0Init,
        	P0=[P0Init[1],P0Init[2],1];
                );

	f=polhomogenise(PolCurveInit,'z);
	if(type(P0)=="t_INT",print("WRONG P0"),	
		if(CheckPointCurve(CCurve,P0)==0,print("\n","BECAREFUL: P0 is not on the curve","\n"));
		);
	if(type(P1)!="t_INT",if(CheckPointCurve(CCurve,P1)==0,print("\n","BECAREFUL: P1 is not on the curve","\n")));
	if(type(P2)!="t_INT",if(CheckPointCurve(CCurve,P2)==0,print("\n","BECAREFUL: P2 is not on the curve","\n")));

	if(CheckP0(CCurve,P0Init),
		flag=0;
		);
	\\change of variables <-> flag==1
	\\if((flagP0==0)||( (#P0Init==3)&&(P0Init[3]==0) ),
	if(flag==0,
		PolCurve=PolCurveInit;
		\\P0=P0Init;
		change=['x,'y,'z];
		);
	\\if(  ((flagP0==0)||( (#P0Init==3)&&(P0Init[3]==0) )) == 0,
	if(flag,
		print("CHANGE OF VARIABLES");
		\\f=polhomogenise(PolCurveInit,'z); we moved this line above
		\\if(#P0Init==3,P0=P0Init,
		\\	P0=[P0Init[1],P0Init[2],1];
		\\	);
	
		b=P0[3];	
		if(P0[3]!=0,
			\\goal: x -> 1/(x-xP0)

			\\x -> x-xP0
			P0=P0/P0[3];
			\\prechange=['x-P0[1]*'z,'y,'z];
			a=P0[1];
			f=subst(f,x,x+P0[1]*z);
			P0[1]=0;

			\\x->1/x
			\\change=['z,'y,'x];
			f=substvec(f,[x,y,z],[z,y,x]);
			P0=[P0[3],P0[2],P0[1]];

			change=['z,'y,'x-a*'z];
			);

		if(b==0,\\then x=0 because P0 is a pole of x and P0=[0,1,0]
			\\goal: x -> y
			change=['y,'x,'z];
			f=substvec(f,[x,y,z],[y,x,z]);
			P0=[P0[2],P0[1],P0[3]];
			);
/*
		if(P0[1..2]!=[0,0],prechange=['x,'y,'z],
			prechange=['x+'z,'y,'z];
                        f=subst(f,x,x-z);
                        P0[1]=1;
                        );
		if(P0[1]!=0,
			change=['x,'y,'z-(1/P0[1])*'x];
			f=subst(f,z,z+(1/P0[1])*x);
			P0[3]=0;
			);
		if(P0[1]==0,
			change=['x,'y,'z-(1/P0[2])*'y];
                        f=subst(f,z,z+(1/P0[2])*y);
                        P0[3]=0;
                        );
*/
		if(P0[1]!=0,P0=P0/P0[1],
			P0=P0/P0[2];
			);
		\\change=substvec(change,[x,y,z],prechange);
		PolCurve=subst(f,z,1);
		
		if(P1Init!=0,
			P1=substvec(change,[x,y,z],P1);
			if(substvec(polhomogenise(PolCurve,'z),[x,y,z],P1)!=0,print("\n","ERROR: after change of variables, new P1 is not right \n"));
			);
		if(P2Init!=0,
                        P2=substvec(change,[x,y,z],P2);
                        if(substvec(polhomogenise(PolCurve,'z),[x,y,z],P2)!=0,print("\n","ERROR: after change of variables, new P2 is not right \n"));
                        );
		print("x -> ",change[1],"  y-> ",change[2],"  z -> ",change[3]);
		);
	\\print("Equation of the curve: ",PolCurve," = 0");
	\\print("x -> ",change[1],"  y-> ",change[2],"  z -> ",change[3]);
	print("P0 = ",P0);
	print("P1 = ",P1);
	print("P2 = ",P2);
	return([PolCurve,P0,change,P1,P2]);
}

InitOmegaArchi(PolCurve)=
\\return a basis of polynomial representing a basis of Omega over Q
{
	my(C,diff,w,Omega);
	C=crvinit(PolCurve,'t,'a);

	\\Check: C is smooth
	if( CheckCurveSingular(C) ,print("\n","WARNING: CURVE NOT SMOOTH","\n"););

	\\INFO CURVE
	crvprint(C);

	diff=C[7][1];
	if(#diff==0,
		printf("\n");
		print("ERROR: GENUS=0");
		printf("\n"););
	w=diff[1];
	Omega=vector(#diff,i,0);
	Omega[1]=1;
	for(k=2,#diff,
		Omega[k]=diff[k]/w);
	print("basis of the regular differentials: ",Omega);
	return([Omega,w,C]);		
}

InitCurveConst(PolCurve)=
{
	my(Omega,w,C,g,d0Initial);
	[Omega,w,C]=InitOmegaArchi(PolCurve);
	\\print("\n","WARNING: OMEGA IS PROBABLY COMPUTED TWICE","\n");
	g=#Omega;
	\\print("GENUS = ",g);
	return([Omega,w,C]);
	\\WARNING: OMEGA IS BASED ON THE MODEL: WE KEEP TRACK OF THE ONE IT HAS BEEN COMPUTED ON
}

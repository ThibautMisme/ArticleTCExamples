\\ \r AaMakdisiMachinery.gp
\\ \r AaTechnique.gp

\\Starting from now, we assume we found a symplectic basis of J[2] over Fq
\\We represent the elements of J[2] as vector of coordinates (Mod 2) in this basis (ei,fi)
\\remark: <X,ei>Weil = Coord(X,fi)

\\B0data=[g,WeilMatrix,vB0,ArfB0]

\\print("ATTENTION: CheckOTC n'est pas optimal"); Ca change rien
\\print("listeOdd: Comment produire tout l'espace sans point random ?"); Ca n'a pas l'air d'être un problème
\\print("CHECKER ARF EN COMPARANT AVEC LES VALEURS DE QUADB0");
\\print("si alpha marche pas, comment recommencer en sauvant de la data ??");
\\print("ATTENTION: si E n'est PAS effectif, il se peut que H°(E)==0 sur Q mais pas sur Fq");

install("PicTors_FrobGen","GGGG")
install("PicTorsSpaceFrobEval","GGGUGG")
install("AllPols","GGUUG")
\\install("BeautyPol","G")
\\install("AlphaNicolas","UGGG")

InitFrob(JFF,g,SBasisJ2FF,WeilMatrix)=
\\InitFrob(JFF,g,SBasisJ2FF,Weildata,WeilMatrix,CCurve)=
\\Frob is the matrix of the Frob over J[2](Fq)
{
	my(Frob,Bk,GB,GBWeil);
	Frob=matrix(2*g,2*g,i,j,Mod(0,2));
	for(j=1,#SBasisJ2FF,
		Bk=SBasisJ2FF[j];
		GB=FrobMakdisi(JFF,Bk);
		\\GBWeil=vector(2*g,i, WeilPairingFF(JFF,CCurve,GB,SBasisJ2FF[i],Weildata))~;
		GBWeil=vector(2*g,i, WeilPairingFFNicolas(JFF,GB,SBasisJ2FF[i]))~;
		for(i=1,#SBasisJ2FF,
			if(i<=g,Frob[i,j]=GBWeil[i+g],
				Frob[i,j]=GBWeil[i-g]));
		);
	return(Frob);
}
/*
InterpolationArf(JFF,ch,WK0tFF,char)=
{
	my(cg,gCoef,g,c,cCoef,gc,h);
	print("h est il dans Fp ?");
	cg=WK0tFF[,1];
	gCoef=matinverseimage(getWn(JFF,2),cg);
	g=JFF[4][2]*gCoef;\\has K0 as zeros at least
	c=ColPointWiseFF(JFF,cg,ch);\\rpz f in H°(2D0)
	cCoef=matinverseimage(getWn(JFF,2),c);
	gc=JFF[4][2]*cCoef;
	h=gc/g;
	return(h);
}


ArfTentative(JFF,W,WK0tFF,P0,char,CCurvep)=
{
	my(d0,g,dim,Wp,Kp,ch,h,D,P,d,RR,Arf,ligneSup);
	d0=getd0(JFF);
	g=getg(JFF);
	dim=2*d0 -g +1;
	Wp=DivAddTM(JFF,W,W,dim);
	dim=d0 -g +1;
	Wp=DivSubTM(JFF,Wp,WK0tFF,dim,getKn(JFF,2));
	Kp=WtoKFF(JFF,Wp);
	ch=DivSubTM(JFF,Wp,getW0(JFF),1,getKn(JFF,2));
	print("dim ch = ",#ch);
	\\ch=Wp[,1] ??
	if(#ch!=1,
		print("WARNING:  Arf, pb)"););
	ch=ch[,1];
	if(char=="Fp",
		h=InterpolationArf(JFF,ch,WK0tFF,char);
		D=FnDiv(CCurvep,h,0);
		P=P0;
		if(#P==2,
			P=[P0[1],P0[2],1];
			);
		d=d0-(g-1);
		ligneSup=[];
		for(k=1,#D[,1],
			if(D[k,2]>0,
				if((D[k,2]%2)!=0,
					print("WARNING: ARF: pn"););
				D[k,2]=D[k,2]/2;
				);
			if(D[k,2]<0,
				ligneSup=concat(ligneSup,[k]);
				);
			if(D[k,1]==P,
				D[k,2]=D[k,2]-d;
				);
			);
		while(#ligneSup!=0,
			k=ligneSup[#ligneSup];
			D=concat(D[1..(k-1),]~,D[(k+1)..(#D[,1]),]~)~;
			);
		print("Arf : E_TC = ",D);
		RR=RiemannRoch(CCurvep,D);
		Arf=Mod(#RR,2);
		);
	return(Arf);
}
*/
/*
RandomPolWTC(JFF,WTC,W2,Vpol,char)=
{
	my(cf,fCoef,fpol);
	cf=WTC*(vector(#WTC,i,random(10))~);
        fCoef=matinverseimage(reduceFFMatrix(JFF,W2),reduceFFCol(JFF,cf));
        fpol=Vpol*fCoef;
	return(fpol);
}

FnDivZeroTM(CCurvep,pol)=
{
	my(D,lignePoles,index);
	D=FnDivTM(CCurvep,pol);
        lignePoles=[];
        for(k=1,#D[,1],
                if(D[k,2]<0,
                        lignePoles=concat(lignePoles,[k]);
                        ););
        while(#lignePoles!=0,
                        index=lignePoles[#lignePoles];
                        D=concat(D[1..(index-1),]~,D[(index+1)..(#D[,1]),]~)~;
			lignePoles=lignePoles[1..(#lignePoles-1)];
                        );
	return(D);
}

IntersectionDiv(D1,D2)=
{
	my(listeindex,mul,D3);
	listeindex=[];
	for(i=1,#D1[,1],
		for(j=1,#D2[,1],
			if(CheckEgalitePointCoord(D1[i,1],D2[j,1]),
				mul=min(D1[i,2],D2[j,2]);
				listeindex=concat(listeindex,[[i,mul]]);
				);););
	D3=matrix(#listeindex,2,i,j,0);
	for(i=1,#listeindex,
		D3[i,1]=D1[listeindex[i][1],1];
		D3[i,2]=listeindex[i][2];
		);
	return(D3);
}

degDiv(D)=
{
	my(deg);
	deg=0;
	for(k=1,#D,
		deg=deg+D[k,2];
		);
	return(deg);
}


ArfBourrin(JFF,WTC,CCurvep,P0,listeBranchP,char)=
{
	my(W2,Vpol,d0,BranchP0,maxdeg,mindeg,minindex,degc,degrandom,g,cg,cCoef,cpol,Dc,fpol,Df,DTC,Pole,Arf);
	W2=getWn(JFF,2);
	Vpol=JFF[4][2];
	d0=getd0(JFF);
	g=getg(JFF);
	BranchP0=listeBranchP[1];

	cg=WTC[,1];
	cCoef=matinverseimage(reduceFFMatrix(JFF,W2),reduceFFCol(JFF,cg));print("cCoef = ",liftFFCol(cCoef));
	cpol=Vpol*cCoef;print("cpol = ",lift(lift(cpol)) );
	Dc=FnDivZeroTM(CCurvep,cpol);print("Dc = ",Dc);
        print("Dc = ",Dc);
	
	fpol=RandomPolWTC(JFF,WTC,W2,Vpol,char);
        Df=FnDivZeroTM(CCurve,fpol);
	degrandom=degP0Pol(fpol,BranchP0);
	if(degrandom==(2*d0),
		degE=d0;
		);
	if(degrandom<(2*d0),
		BranchP0=listeBranchP[1];
        	maxdeg=0;
        	mindeg=2*d0;
        	minindex=1;
        	for(k=1,#WTC,
                	cg=WTC[,k];
                	cCoef=matinverseimage(reduceFFMatrix(JFF,W2),reduceFFCol(JFF,cg));
                	cpol=Vpol*cCoef;
                	degc=-degP0Pol(cpol,BranchP0);print("degc = ",degc);
                	maxdeg=max(maxdeg,degc);
                	if(degc<mindeg,
                        	mindeg=degc;
                        	minindex=k;
                        	);
                	);
        	degE=d0-(2*d0-maxdeg);
		print("maxdeg = ",maxdeg," mindeg = ",mindeg," minindex = ",minindex," degE = ",degE," d0 = ",getd0(JFF));
		);

	DTC=IntersectionDiv(Dc,Df);print("DTC = ",DTC);
	while(degDiv(DTC)>degE,
		fpol=RandomPolWTC(JFF,WTC,W2,Vpol,char);
        	Df=FnDivZeroTM(CCurve,fpol);
        	DTC=IntersectionDiv(DTC,Df);print("DTC = ",DTC);
		);
	if(#DTC!=d0,
		print("ArfBourrin: pb");
		);
	
	Pole=matrix(1,2,i,j,0);
	Pole[1,1]=P0;
	if(#P0==2,
		Pole[1,1]=[P0[1],P0[2],1];
		);
	Pole[1,2]=-(d0-(g-1));
	DTC=concat(DTC~,Pole~)~;print("final DTC = ",DTC);
	
	Arf=#RiemannRochTM(CCurvep,DTC);print("dim H°(B) = ",Arf);
	Arf=Mod(Arf,2);
	return(Arf);
}

Arf2(JFF,WTC,WE1,WArf)=
\\Warf=H°( (d0-g+1)P0 )
{
	my(d0,g,V,S1p,K1p,S1,s1,sV,KsV,S2p,K2p,W,Arf);
	d0=getd0(JFF);
	g=getg(JFF);
	V=getWn(JFF,2);
	S1p=DivAddTM(JFF,WTC,WE1,(4*d0 -d0 -(d0-g)) -g +1);
	K1p=WtoKFF(JFF,S1p);
	S1=DivSubTM(JFF,S1p,V,1,K1p);\\generically of dim 1. And we only need one function of this space anyway
	s1=S1[,1];
	sV=MulPointWiseFF(JFF,s1,V);
	KsV=WtoKFF(JFF,sV);
	S2p=DivSubTM(JFF,sV,WTC,d0-g+1,KsV);\\H°(2D0 -E1 -Ex)~H°(E) si WTC=H°(2D0 -E)
	K2p=WtoKFF(JFF,S2p);
	W=DivSubDimTM(JFF,S2p,Warf,K2p);
	Arf=Mod(#W,2);
	return(Arf);
}
*/

/*
Arf3(JFF,WTC,WArf,BranchP0)=
\\WArf=H°( (d0-g+1)P0 )
{
	my(s,sV,KsV,WEs,KEs,W,Arf);
	s=WTC[,1];print("deg P0 s = ",-degP0Col(JFF,s,BranchP0));
	while((-degP0Col(JFF,s,BranchP0))<(2*getd0(JFF)),
		s=WTC*vector(#WTC,i,random(10))~;print("deg P0 s = ",-degP0Col(JFF,s,BranchP0));
		);

	sV=MulPointWiseFF(JFF,s,getWn(JFF,2));\\H°(2D0 - (D+Es -2D0))=H°(4D0-D-Es)
	\\KsV=WtoKFF(JFF,sV); NON
	WEs=DivSubTM(JFF,sV,WTC,getd0(JFF)-getg(JFF)+1,getKn(JFF,2));\\H°(2D0 -Es) avec 2D0-Es~D tq WTC=H°(2D0-D)
	KEs=WtoKFF(JFF,WEs);

	W=DivSubDimKFF(JFF,WEs,WArf,KEs);\\H°(2D0-Es -dP0) ~ H°(D-dP0)
	\\ATTENTION: H°(dP0) n'est peut-etre pas base-point free
	\\W=H°(2D0-Es -dP0) = {v in WEs| deg_P0(v) <= 2d0-d} IF Supp(Es) et P0 disjoints

	Arf=Mod(#W,2);print("dim H°(Theta Caracteristic) = ",#W," W = ",W);
	return(Arf);
}
*/

Arf4(JFF,WTC,WArf,CCurve,P0)=
\\WArf=H°(2D0 -dP0) = H°(D0 +(g-1)P0)
\\d=d0 -g +1
\\because d+(g-1)=d0
{
        my(s,sV,KsV,WEs,KEs,Wp,W,Arf);
	if(default(debug)>1,print("Arf 4"));
        s=WTC[,1];\\print("deg P0 s = ",-degP0Col(JFF,s,CCurve,P0));
        while((-degP0Col(JFF,s,CCurve,P0))<(2*getd0(JFF)),
                s=WTC*vector(#WTC,i,random(10))~;\\print("deg P0 s = ",-degP0Col(JFF,s,CCurve,P0));
                );

        sV=MulPointWiseFF(JFF,s,getWn(JFF,2));\\H°(2D0 - (D+Es -2D0))=H°(4D0-D-Es)
        \\KsV=WtoKFF(JFF,sV); NON
        WEs=DivSubTM(JFF,sV,WTC,getd0(JFF)-getg(JFF)+1,getKn(JFF,2));\\H°(2D0 -Es) avec 2D0-Es~D tq WTC=H°(2D0-D)
        KEs=WtoKFF(JFF,WEs);

	Wp=DivAddTM(JFF,WEs,WArf,2*getd0(JFF));
        W=DivSubDimKFF(JFF,Wp,getWn(JFF,2),KEs);\\H°(2D0-Es -dP0) ~ H°(D-dP0)
	\\print(W); 
        Arf=Mod(#W,2);if(default(debug)>1,print("dim H°(Theta Caracteristic) = ",#W));
        return(Arf);
}

Arf(JFF,WTC,Arfdata)={return(Arf4(JFF,WTC,Arfdata[1],Arfdata[2],Arfdata[3]));}

/*
QuadCanonique(g,v)=
\\g = genus = #(ei) = #fi
{
	my(res);
	res=0;
	for(k=1,g,
		res=res + v[k]*v[k+g];);
	return(res);
}
*/

/*
InitQuadB0(JFF,SBasisJ2FF,WB0tFF)=
\\QuadB0=QuadCanonique + <,vB0>Weil
{
	my(g,ArfB0,listeValQuadB0,vB0);
	g=getg(JFF);
	ArfB0=Arf(JFF,WB0tFF,getW0(JFF));
	listeValQuadB0=vector(#SBasisJ2FF,i,Arf(JFF,WB0tFF,SBasisJ2FF[i])+ArfB0 );
	vB0=vector(2*g,i,0)~;
	for(k=1,g,
		vB0[k]=listeValQuadB0[g+k] - Mod(0,2) ;);\\QuadCanonique is zero over a symplectic basis
	for(k=1,g,
		vB0[g+k]=listeValQuadB0[k] - Mod(0,2) ;);
	return([vB0,ArfB0]);
}*/



QuadB0Old(B0data,v)=
\\WeilMatrix is supposed to be the symplectic Matrix
{
	\\my(g,WeilMatrix,vB0,ArfB0);
	\\[g,WeilMatrix,vB0,ArfB0]=B0data;
	return(QuadCanonique(B0data[1],v) + mattranspose(B0data[3])*B0data[2]*v);
}



CheckOTCOld(B0data,v)=
{
	my(g,WeilMatrix,vB0,ArfB0,res);
	[g,WeilMatrix,vB0,ArfB0]=B0data;
	if(lift(ArfB0)==1,\\We look for the zeros of QuadB0
                res=(lift(QuadB0Old(B0data,v))==0);
                );
        if(lift(ArfB0)==0,\\We look for the NON-zeros of QuadB0
                res=(lift(QuadB0Old(B0data,v))==1);
                );
	return(res);
}


CheckOTC(B0data,v)=
\\1 si OTC
{
	my(res);
	\\[g,WeilMatrix,vB0,ArfB0]=B0data;
	\\res=-1;
	\\What we want to do:
/*
	if(lift(ArfB0)==1,\\We look for the zeros of QuadB0
		res=(lift(QuadB0(B0data,v))==0);
		);
	if(lift(ArfB0)==0,\\We look for the NON-zeros of QuadB0
		res=(lift(QuadB0(B0data,v))==1);
		);
*/
	\\It is equivalent to:
	res = (lift(QuadB0(B0data,v)))!= (lift(B0data[4])) ;
	\\print("A EFFACER A TERME: CheckOTC = CheckOTCOLD : ",res == CheckOTCOld(B0data,v));
	return(res);
}

InitMatA(Frob)={return(Mod(matid(#Frob),2)-Frob);}

InitSymplecticMatrix(g)=
{
	my(Id,Zero,WMat);
	Id=Mod(matid(g),2);
	Zero=matrix(g,g,i,j,Mod(0,2));
	WMat=concat(Zero,Id);
	WMat=concat(WMat~,concat(Id,Zero)~)~;
	return(WMat);
}

/*
IntBasis2(ninit,g)=
{
	my(m,n,decomp,eucl);
	m=2*g;
	n=ninit;
	decomp=vector(m,i,0);
	for(k=1,m,
		eucl=n%2;
		decomp[k]=Mod(eucl,2);
		n=(n-eucl)/2;
		);
	return(decomp);
}

InitlisteOdd(g,B0data)=
\\Coef of Point P of J[2] in the basis SBasisJ2FF s.t. B0 + P are all the OTC
{
	my(CardOTC,listeOdd,k,otc,v);
	CardOTC=(2^(g-1))*(2^g-1);
	listeOdd=vector(CardOTC,i,0);
	k=0;
	otc=0;
	while((otc<CardOTC)&&(k<2*g),
		k=k+1;
		v=IntBasis2(k,g)~;
		if(CheckOTC(B0data,v),
			otc=otc+1;
			listeOdd[otc]=v;
			);
		);
	if(otc<CardOTC,
		print("WARNING: ALL OTC NOT FOUND"););
	return(listeOdd);
}
*/

RecoverlistOTCNt(listeOrbiteOdd,g)=
{
	my(CardOTC,listOTCNt,index);
        CardOTC=(2^(g-1))*(2^g-1);
	listOTCNt=vector(CardOTC,i,0);
	index=0;
	for(k=1,#listeOrbiteOdd,
		for(l=1,#listeOrbiteOdd[k],
			index=index+1;
			listOTCNt[index]=DicoTCFtoNt(listeOrbiteOdd[k][l]);
			);
		);
	listOTCNt=Set(listOTCNt);
	return(listOTCNt);		
}


InitlisteOrbiteOddCore(B0data,Frob,flagTC)=
\\return the list of the orbits of J[2] corresponding to OTC (flagTC=0) or to ETC (flagTC=1) under the Frob
{
	my(g,listeOrbiteOdd,k,CardOTC,listOTCNt,counter,v,Orbitev,Galv);
	print("\n");
	if(flagTC,print("	EVEN TC computation"),
		print("		ODD TC computation"));	
	listeOrbiteOdd=[];
	k=0;
	g=B0data[1];
	if(flagTC,CardOTC=(2^(g-1))*(2^g+1),\\cardOTC is the cardinal of ETC in that case
		CardOTC=(2^(g-1))*(2^g-1));
	if(flagTC,print("CardETC = ",CardOTC),
		print("CardOTC = ",CardOTC));
	listOTCNt=vector(CardOTC,i,0);
	while(k<CardOTC,

		print1("new orbit: #",k,"   ");
		counter=0;

		v=vector(2*g,i,Mod(random(2),2))~;
		\\while((CheckOTCOld(B0data,v)==0)||IsInListNested(listeOrbiteOdd,v),
		while((CheckOTCOld(B0data,v)==flagTC)||IsInListNested(listeOrbiteOdd,v),
			v=vector(2*g,i,Mod(random(2),2))~;

			counter=counter +1;

			);

		if(default(debug)>2,print("counter random point : ",counter));

		Orbitev=[v];
		\\kold=k;
		k=k+1;
		listOTCNt[k]=DicoTCFtoNt(v);
		Galv=Frob*v;
		while(Galv!=v,
			Orbitev=concat(Orbitev,[Galv]);
			k=k+1;
			listOTCNt[k]=DicoTCFtoNt(Galv);
			Galv=Frob*Galv;
			);
		\\print("check k: ",k==(kold+#Orbitev));
		print("size of the added Fq-orbit: ",#Orbitev);
		listeOrbiteOdd=concat(listeOrbiteOdd,[Orbitev]);
		);
	if(k!=CardOTC,print("ERROR InitlisteOrbiteOdd: Not all (Odd or Even) TC have been found"));
	listOTCNt=Set(listOTCNt);
	return([listeOrbiteOdd,listOTCNt]);
}

InitlisteOrbiteOdd(B0data,Frob,flagTC)=
\\return the list of the orbits of J[2] corresponding to OTC under the Frob
\\if flagTC=1, we return as well the list of orbits of ETC. Then listeOrbiteOdd=[ listeOrbiteOdd , listeOrbiteEven ]
{
	my(listeOrbiteOdd,listOTCNt,listeOrbiteEven,listETCNt);
	[listeOrbiteOdd,listOTCNt]=InitlisteOrbiteOddCore(B0data,Frob,0);
	if(flagTC,
		[listeOrbiteEven,listETCNt]=InitlisteOrbiteOddCore(B0data,Frob,1);
		listeOrbiteOdd=[listeOrbiteOdd,listeOrbiteEven];
		listOTCNt=[listOTCNt,listETCNt];
		);
	return([listeOrbiteOdd,listOTCNt]);
}
		
InitlisteOrbiteJ2(g,Frob)=
{
	my(listeOrbiteJ2,k,CardJ2,v,Orbitev,Galv);
	print("size J[2] = ",2^(2*g));
	listeOrbiteJ2=[];
	k=0;
	CardJ2=2^(2*g);
	while(k<CardJ2,
		v=vector(2*g,i,Mod(random(2),2))~;
		while(IsInListNested(listeOrbiteJ2,v),
                        v=vector(2*g,i,Mod(random(2),2))~;
                        );
		Orbitev=[v];
                Galv=Frob*v;
                while(Galv!=v,
                        Orbitev=concat(Orbitev,[Galv]);
                        Galv=Frob*Galv;
                        );
                k=k+#Orbitev;
                print("size listeOrbiteJ2 = ",k," size orbit over Fq = ",#Orbitev);
		listeOrbiteJ2=concat(listeOrbiteJ2,[Orbitev]);
		);
	return(listeOrbiteJ2);
}
		



SuppImA(M,A)=
{
	my(dim,Jp,Jp2,S);
	dim=#M[,1];
	Jp=matker(A);
	Jp2=matker(A*A);
	if(#Jp2==#Jp,\\then J[2](Fp^2) == J[2](Fp) and J[2](Fp) is a supplement of ImA 
		S=Jp;
		print("J[2](Fp^2) == J[2](Fp)    dim =",#Jp);
		if(matrank(concat(M,S))!=dim,
			print("WARNING: S IS NOT A SUPPLEMENT OF ImA");
			S=matsupplement(M)[,(#M+1)..dim];
			);
		);
	if(#Jp2!=#Jp,
		S=matsupplement(M)[,(#M+1)..dim];
		\\V=matsupplement(A);
		\\S=V[,(#A+1)..#V];
		);
	return(S);
}

Checkv0(testArf,B0data)=
{
	my(flag);
	flag=1;
	for(k=1,#testArf,
		if( QuadB0(B0data,testArf[k][2]) != (testArf[k][1] + B0data[4]), \\we check that q_B0(B-B0)== arf + ArfB0 indeed (for every OTC such that we have computed Arf)
			print("ERROR: Arf computation");
			flag=0;
			); 
		);
	return(flag);
}

Initv0(g,A,M,S,Frob,JFF,WB0tFF,SBasisJ2FF,WeilMatrix,Arfdata,ArfB0)=
{
	my(Weilv0,m,v,c,testArf,bCoef,index,b,WTC,Bt,v0);
	Weilv0=vector(2*g,i,0)~;\\the vector <v0,v> when v ranges through a basis of M+S
	m=#M;
	for(k=1,m,
		v=M[,k];
		c=matinverseimage(A,v);
		Weilv0[k]=QuadCanonique(g,Frob*c)-QuadCanonique(g,c);
		);
	/*
	if(#S!=0,
		fpol=Mod(JFF[1],getpMak(JFF));
        	CCurvep=CrvInit(fpol);
		);
	*/
	testArf=[];
	print("Number of components of v0 to compute via Direct Arf Computation: ",#S);
	for(k=1,#S,
		bCoef=S[,k];
		index=FirstNonZero(lift(bCoef));
		b=SBasisJ2FF[index];
		for(i=index+1,#SBasisJ2FF,
			if(lift(bCoef[i])!=0,
				b=SumJ(JFF,b,SBasisJ2FF[i]);
				);
			);
		WTC=SumJ(JFF,WB0tFF,b);
		arf=Arf(JFF,WTC,Arfdata);
		Weilv0[m+k]=arf-QuadCanonique(g,bCoef)-ArfB0;
		testArf=concat(testArf,[[arf,bCoef]]);
		);
	Bt=concat(M,S);
	v0=matsolve(mattranspose(Bt)*WeilMatrix,Weilv0);
	return([v0,testArf]);
}	











\\************************************************************** LIFT BASIS ******************************************
\\********************************************************************************************************************

LiftBasisJ2FF(J,BasisJ2FF)=
{
\\return(vector(#BasisJ2FF,i,LiftTorsion(J,BasisJ2FF[i])));
	my(BasisJ2GAdic);
        BasisJ2GAdic=vector(#BasisJ2FF,i,0);
        for(i=1,#BasisJ2GAdic,
                BasisJ2GAdic[i]=LiftTorsion(J,BasisJ2FF[i]);
                print("Lift Basis: ",i,"/",2*getg(J));
                );
        return(BasisJ2GAdic);
}







export(LiftTorsion)

LiftBasisJ2F2FrobGenLiftLoop(J,JFF,BasisJ2FF,MFrob)=
{
	my(PicTorsFrobGendata,Gen,CG,MFrobcan,GenGAdic);
	PicTorsFrobGendata=PicTors_FrobGen(JFF,2,BasisJ2FF,MFrob);
        [Gen,CG,MFrobcan]=PicTorsFrobGendata;
        GenGAdic=vector(#Gen,k,0);
	print("\n");
        print("Lift of the F2[Frob]-generators of J[2]:");

/*
        for(k=1,#Gen,
                GenGAdic[k]=LiftTorsion(J,Gen[k]);
                print("         ",k,"/",#Gen" done"); 
                );
	A=GenGAdic[1];
*/
	parfor(k=1,#Gen,
		print("		",k,"/",#Gen," beginning");
		LiftTorsion(J,Gen[k]),
			r,GenGAdic[k]=r);

	\\B=GenGAdic[1];print("A == B ??",CheckEqual(J,A,B));

        PicTorsFrobGendata[1]=GenGAdic;
	PicTorsFrobGendata[3]=MFrob;
	return(PicTorsFrobGendata);
}

RecoverBasisGAdic(J,PicTorsFrobGendata,dimBasisJ2,MFrob)=
{
	my(GenGAdic,CG,MFrobcan,M,GenLinGAdic,flag,dim,vgen,v,Wv,BasisJ2GAdic,coef);
	[GenGAdic,CG,MFrobcan]=PicTorsFrobGendata;
	M=matrix(#CG[1],#CG,i,j,Mod(CG[j][i],2));
        GenLinGAdic=GenGAdic;\\will become a linear spanning family
        flag=1;
        dim=dimBasisJ2;
        for(k=1,#GenGAdic,
		if(flag,
                	vgen=M[,k];
                	v=MFrob*vgen;
                	Wv=picfrob(J,GenGAdic[k]);
                	while(v!=vgen,
                        	if(flag,
                                	M=concat(M,v);
                                	GenLinGAdic=concat(GenLinGAdic,[Wv]);
                                	flag=(matrank(M)!=dim);\\print("k = ",k," flag = ",flag);
                                	Wv=picfrob(J,Wv);
                                	);
                        	v=MFrob*v;
                        	);
                	);
		);
        if(flag,print("\n error LiftBasisJ2F2FrobGen: we couldn't get a linear spanning family. dim =",dim," dim family =",matrank(M),"\n"));
        BasisJ2GAdic=vector(dim,i,0);
        for(k=1,dim,
                coef=matinverseimage(M,vector(dim,i,Mod(i==k,2))~);
                BasisJ2GAdic[k]=piclc(J,coef,GenLinGAdic);
                );
	return(BasisJ2GAdic);
}


LiftBasisJ2F2FrobGen(J,JFF,BasisJ2FF,MFrob)=
{
	my(PicTorsFrobGendata,BasisJ2GAdic);
	PicTorsFrobGendata=LiftBasisJ2F2FrobGenLiftLoop(J,JFF,BasisJ2FF,MFrob);
	BasisJ2GAdic=RecoverBasisGAdic(J,PicTorsFrobGendata,#BasisJ2FF,MFrob);
	return([BasisJ2GAdic,PicTorsFrobGendata]);
}

\\*****************************************************************************************************








InitWArf(JFF,CCurvep,P0)=
{
	return(PoltoMakdisi(JFF,RiemannRochTM(CCurvep,[P0,getd0(JFF) +getg(JFF) -1])));
}

/*
LiftTC2(J,WTCinit,WK0tGAdic,CardJq,CardJp)=
{
	my(p,g,logemax,WTC,e,emax,pe,CardJqe,CardJpe,Mt,WK0,u,v,d,N,U,Jinter,WTCold);
	print("WARNING: LiftTC n'est pas assuré de marcher");
	p=getpMak(J);
	\\print("ON PROFITE DU FAIT QUE WK0tGAdic ET WB0tGAdic SONT DES POINTS DE J(Fp) POUR DIMINUER U ET N EN LES REDUISANT MOD CardJpe");
	g=getg(J);
	emax=geteJ(J);
	print("LiftTC: beginning. emax = ",emax);
	logemax=logint(emax,2);
	WTC=WTCinit;
	for(k=1,logemax,
		e=2^k;\\precision we aim for
		Jinter=RedJ(J,e);
		print("LiftTC: e = ",e);
		pe=p^e;
		print("pe = ",pe);
		CardJqe=CardJq*(p^(g*(e-1)));print("CardJqe = ",CardJqe);
		CardJpe=CardJp*(p^(g*(e-1)));print("CardJpe = ",CardJpe);
		WTCold=WTC;
		Mt=LiftRandom(Jinter,2^(k-1),WTC);\\a point of J(Zq/p^e) to the precision e such that it is equal to the matrix WTC to the precision e/2
		WK0=RedPointJ(J,WK0tGAdic,e);
		\\[u,v,d]=gcdext(2,pe);
		[u,v,d]=gcdext(pe,CardJp);
		\\N=(v*pe)%CardJp;print("N = ",N);
		if(abs(u*pe)>=(u*pe)%CardJqe,
			N=(u*pe)%CardJqe;
			);
		if(abs(u*pe)<(u*pe)%CardJqe,
			N=u*pe;
			);
		print("N = ",N);print("u*pe = ",u*pe);print("(u*pe)%CardJqe = ",(u*pe)%CardJqe);
		\\U=u%CardJp;print("U = ",U);
		\\WTC=SumJ(Jinter,MulInt(Jinter,Mt,N),MulInt(Jinter,WK0,U));\\the good lift
		WTC=MulInt(Jinter,Mt,N);
		print("Test WB0tGAdic : ",TestWB0tGAdicK0(Jinter,WTC,WK0)," et ",TestWB0tGAdicRed(Jinter,WTC,WTCold,2^(k-1)));
		);
	if(emax!=e,
		e=emax;print("LiftTC: e = ",e);
                pe=p^e;print("pe = ",pe);
                CardJqe=CardJq*(p^(g*(e-1)));print("cardJqe = ",CardJqe);
		WTCold=WTC;
                Mt=LiftRandom(J,2^logemax,WTC);\\a point of J(Zq/p^e) to the precision e such that it is equal to the matrix WTC to the precision e/2
                WK0=WK0tGAdic;
                [u,v,d]=gcdext(pe,CardJp);
                if(abs(u*pe)>=(u*pe)%CardJqe,
                        N=(u*pe)%CardJqe;
                        );
                if(abs(u*pe)<(u*pe)%CardJqe,
                        N=u*pe;
                        );
		print("N = ",N);print("u*pe = ",u*pe);print("(u*pe)%CardJqe = ",(u*pe)%CardJqe);
                WTC=MulInt(J,Mt,N);
		print("Test WB0tGAdic : ",TestWB0tGAdicK0(J,WTC,WK0)," et ",TestWB0tGAdicRed(J,WTC,WTCold,2^logemax));
		);
	return(WTC);
}
*/

LiftTC1method(J,WTCinit,WK0tGAdic,CardJq,CardJp,emax)=
{
        my(p,g,logemax,WTC,e,pe,CardJqe,CardJpe,Mt,WK0,u,v,d,N,U,Jinter,WTCold);
        p=getpMak(J);
        \\print("ON PROFITE DU FAIT QUE WK0tGAdic ET WB0tGAdic SONT DES POINTS DE J(Fp) POUR DIMINUER U ET N EN LES REDUISANT MOD CardJpe");
        g=getg(J);
        \\emax=geteJ(J);
        \\print("\n","LiftTC: beginning");
	if(default(debug)>1,print("first method"));
        logemax=logint(emax,2);
        WTC=WTCinit;
        for(k=1,logemax,
                if(default(debug)>1,print("\n"));
		e=2^k;\\precision we aim for
                Jinter=RedJ(J,e);
                if(default(debug)>1,print("LiftTC: e = ",e));
                pe=p^e;
                if(default(debug)>1,print("pe = ",pe));
                CardJqe=CardJq*(p^(g*(e-1)));if(default(debug)>1,print("CardJqe = ",CardJqe));
                CardJpe=CardJp*(p^(g*(e-1)));if(default(debug)>1,print("CardJpe = ",CardJpe));
                WTCold=WTC;
                Mt=LiftRandom(Jinter,2^(k-1),WTC);\\a point of J(Zq/p^e) to the precision e such that it is equal to the matrix WTC to the precision e/2
                WK0=RedPointJ(J,WK0tGAdic,e);
                [u,v,d]=gcdext(2,pe);
                \\N=(v*pe)%CardJqe;print("N = ",N);
                if(abs(v*pe)>=(v*pe)%CardJqe,
                        N=(v*pe)%CardJqe;
                        );
                if(abs(v*pe)<(v*pe)%CardJqe,
                        N=v*pe;
                        );
		if(default(debug)>1,print("N = ",N);print("v*pe = ",v*pe);print("(v*pe)%CardJqe = ",(v*pe)%CardJqe););
                \\U=u%CardJp;print("U = ",U);
		if(abs(u)>=u%CardJpe,
			U=u%CardJpe;
			);
		if(abs(u)<u%CardJpe,
			U=u;
			);
		if(default(debug)>1,print("U = ",U));
                WTC=SumJ(Jinter,MulInt(Jinter,Mt,N),MulInt(Jinter,WK0,U));\\the good lift
                if(default(debug)>1,print("Test WB0tGAdic : ",TestWB0tGAdicK0(Jinter,WTC,WK0)," and ",TestWB0tGAdicRed(Jinter,WTC,WTCold,2^(k-1))););
                );
        if(emax!=e,
                e=emax;if(default(debug)>1,print("LiftTC: e = ",e));
                pe=p^e;if(default(debug)>1,print("pe = ",pe));
                CardJqe=CardJq*(p^(g*(e-1)));if(default(debug)>1,print("cardJqe = ",CardJqe));
                WTCold=WTC;
                Mt=LiftRandom(J,2^logemax,WTC);\\a point of J(Zq/p^e) to the precision e such that it is equal to the matrix WTC to the precision e/2
                WK0=WK0tGAdic;
                [u,v,d]=gcdext(2,pe);
                \\N=(v*pe)%CardJqe;print("N = ",N);
                if(abs(v*pe)>=(v*pe)%CardJqe,
                        N=(v*pe)%CardJqe;
                        );
                if(abs(v*pe)<(v*pe)%CardJqe,
                        N=v*pe;
                        );
                if(default(debug)>1,print("N = ",N);print("v*pe = ",v*pe);print("(v*pe)%CardJqe = ",(v*pe)%CardJqe););
                \\U=u%CardJp;print("U = ",U);
                if(abs(u)>=u%CardJpe,
                        U=u%CardJpe;
                        );
                if(abs(u)<u%CardJpe,
                        U=u;
                        );
                if(default(debug)>1,print("U = ",U));
                WTC=SumJ(J,MulInt(J,Mt,N),MulInt(J,WK0,U));\\the good lift
                if(default(debug)>1,print("Test WB0tGAdic : ",TestWB0tGAdicK0(J,WTC,WK0)," and ",TestWB0tGAdicRed(J,WTC,WTCold,2^logemax)););
                );
        return(WTC);
}


LiftTC2method(J,WTC,e,WK0tGAdic)=
{
	return(piclift(J,WTC,e,2,0,WK0tGAdic));
}

LiftTC(Jinit,WTCinit,WK0tGAdic,CardJq,CardJp)=
{
	my(J,emiddle,WTCold,WTC);
	J=Jinit;
	if(default(debug)>1,print("\n","LiftTC: beginning. Going toward ",geteJ(J)));
	emiddle=2;
	WTCold=LiftTC1method(J,WTCinit,WK0tGAdic,CardJq,CardJp,emiddle);
	if(default(debug)>1,print("\n","second method"));
	WTC=LiftTC2method(J,WTCold,emiddle,WK0tGAdic);
	
	\\if((flagP0==0)&&((#P0==2)||((#P0==3)&&(P0[3]==1))),
	if(TestWB0tGAdicK0(J,WTC,WK0tGAdic)==0,
		while((geteJ(J)>10)&&(TestWB0tGAdicK0(J,WTC,WK0tGAdic)==0),
			print("loss of precision (-5) because P0 is affine");
                        J[7]=geteJ(J)-5;
                        );
                print("new precision: ",geteJ(J));
                );

	print("Test WB0tGAdic : ",TestWB0tGAdicK0(J,WTC,WK0tGAdic)," and ",TestWB0tGAdicRed(J,WTC,WTCold,emiddle));
	return([J,WTC]);
}


















\\**************************************************************** CHI *******************************************************
\\****************************************************************************************************************************





\\**************************************************************** ALPHA ****************************************************

AlphaMapNicolas(J,W,i,j)=
{
	my(alphaNicolas,alpha,alp);
	alphaNicolas=AlphaMap(J,W);
	alpha=alphaNicolas[1,1];
	alp=divpadic(J,alpha[i],alpha[j]);
	return(alp);
}

AlphaMapTHIBAUT(J,W,i,j)=
{
	my(d0,g,V,KV,WE1,WE2,S1p,S1,s1,s1V,S2p,S2pp,S2,s2,alpha,flag,alp,p);
	d0=getd0(J);
	g=getg(J);
	V=getWn(J,2);
	KV=getKn(J,2);
	WE1=J[14][1][1];
	WE2=J[14][2][1];
	S1p=DivAddTM(J,W,WE1,2*d0+1);
	S1=DivSubTM(J,S1p,V,1,KV);
	s1=S1[,1];print("s1 == [] (0 expected) ? ",s1==[]);
	s1V=mulPointWise(J,s1,V);
	S2p=DivSubTM(J,s1V,W,d0-g+1,KV);
	S2pp=DivAddTM(J,S2p,WE2,2*d0+1);
	S2=DivSubTM(J,S2pp,V,1,KV);
	s2=S2[,1];print("s2 == [] (0 expected) ? ",s2==[]);

	alpha=ColtoCoefGAdicV2(J,s2);\\print("alpha = ",alpha);

	alp=0;	
	flag=1;
	p=getpMak(J);
	if(liftall(Mod(alpha[i],p))==0,
		flag=2;);
	if(liftall(Mod(alpha[j],p))==0,
		flag=3;);
	\\print("flag = ",flag);
	if(flag==1,
		alp=divpadic(J,alpha[i],alpha[j]);
		\\alp=divpadic(J,alpha[i]+alpha[j+1],alpha[j]);
		);
	if(flag!=1,
		\\print("AlphaMap ended up with a 0 : ",alpha);
		);
	\\alp=rpzQqtoJ(rpzJtoQq(J,alpha[i])/rpzJtoQq(J,alpha[j]));
	return([alp,flag]);
}


InitWE(J,CCurve,d0,g,P0,P,flagW)=
{
	my(E,Epol,WE,res);
        if(CheckEgalitePointCoord(P0,P)||(P==0),E=[P0,d0+g],
                E=matrix(2,2,i,j,0);
                E[1,]=[P0,2*d0];
                E[2,]=[P,-(d0-g)];
                );
        Epol=RiemannRochTM(CCurve,E);
	res=Epol;
	if(flagW,
        	WE=PoltoMakdisi(J,Epol);
		res=[Epol,WE];
		);
	return(res);
}

StartingAlpha(J,CCurve,P0,P1,P2)=
{
	my(E1,E1pol,WE1,E2pol,WE2,MU,d0,g);
	d0=getd0(J);
	g=getg(J);
	[E1pol,WE1]=InitWE(J,CCurve,d0,g,P0,P1,1);
	if(P2!=0,[E2pol,WE2]=InitWE(J,CCurve,d0,g,P0,P2,1),
		E2pol=E1pol;
		WE2=WE1;
		);
	MU=0;
	return([E1pol,E2pol,[[WE1],[WE2],MU]]);
}

/*
ChiOne(J,W,CCurve,P0)=
{
	my(flag,res);
	J[4]=[J[4][1],J[4][2],[],[]];
	[J[4][3],J[4][4],J[14]]=StartingAlpha(J,CCurve,P0);
	res=AlphaMapTHIBAUT(J,W);
	
	flag=0;
	while(flag,
		[J[4][3],J[4][4],J[14]]=NewAlpha();
		res=AlphaMap(J,W);
		flag=0;
		);
	return(res);
}*/

\\*********************************************************************************************************************************

OrbiteElt(J,alp)=
{
        my(listealp,alpFrob,k);
        \\alp=rpzQqtoJ(rpzJtoQq(J,alpha[2])/rpzJtoQq(J,alpha[1]));print("poldegree alp = ",poldegree(alp));
        listealp=[alp];\\print("alp = ",alp);
        alpFrob=FrobQq(J,alp);
	k=1;
        while(alpFrob!=alp,
                listealp=concat(listealp,[alpFrob]);
                alpFrob=FrobQq(J,alpFrob);\\print("alpFrob = ",alpFrob);
		k=k+1;
		);
        return(listealp);
}


\\***************************************************************** CHI J[2] *****************************************************************

InitlisteRpzOrbiteJ2(J,listeOrbiteJ2,BasisJ2GAdic)=
{
	my(listeRpzOrbiteJ2,vectornul,Orbite,coef,WCoef,W);
	listeRpzOrbiteJ2=vector(#listeOrbiteJ2,i,0);
	vectornul=vector(#BasisJ2GAdic,i,Mod(0,2));
	for(k=1,#listeOrbiteJ2,
                Orbite=listeOrbiteJ2[k];
                coef=Orbite[1];\\print("coef = ",lift(coef));
		if(coef!=vectornul,
			W=piclc(J,coef,BasisJ2GAdic);
			listeRpzOrbiteJ2[k]=W;
			);
		/*
                WCoef=CoefNonZero(coef);\\print("WCoef = ",WCoef);
                if(WCoef!=[],\\we remove 0
                        W=BasisJ2GAdic[WCoef[1]];
                        for(i=2,#WCoef,
                                W=SumJTorsion(J,W,BasisJ2GAdic[WCoef[i]]);
                                );
			listeRpzOrbiteJ2[k]=W;
			);
		*/
		print("		orbit ",k,"/",#listeOrbiteJ2," done");
		);
	return(listeRpzOrbiteJ2);
}



TentativeChiJ2(J,listeRpzOrbiteJ2,listeOrbiteJ2,iAlpha,jAlpha)=
{
        my(listeOrbiteAlpha,flagglobal,W,alpha,flagalpha,Orbitealpha,F);
	listeOrbiteAlpha=vector(#listeRpzOrbiteJ2,i,0);
	flagglobal=1;
	for(k=1,#listeOrbiteAlpha,
		if(flagglobal==1,
			W=listeRpzOrbiteJ2[k];
			if(W!=0,\\we remove 0
				[alpha,flagalpha]=AlphaMapTHIBAUT(J,W,iAlpha,jAlpha);
				if(flagalpha!=1,
					flagglobal=flagalpha;);
				\\alpha=AlphaMapNicolas(J,W,iAlpha,jAlpha);
				if(flagglobal==1,
					Orbitealpha=OrbiteElt(J,alpha);
					listeOrbiteAlpha[k]=Orbitealpha;
					print("\n");
					print("size Orbit = ",#listeOrbiteJ2[k]);
					print("size Orbit alpha = ",#Orbitealpha);
					);
				);
			);
		);
	F=1;
	if(flagglobal==1,
		for(i=1,#listeOrbiteAlpha,
			for(j=1,#listeOrbiteAlpha[i],
				F=F*(x-rpzJtoQq(J,listeOrbiteAlpha[i][j]));
				););
		);
	return([F,listeOrbiteAlpha,flagglobal]);		
}

CleanChi(Finit)=
{
	my(F,res,coef);
	F=Finit;
	F=bestappr(F);
	if(F==[],
		print("fail Chi: ",F);
		res=0;
		);
	if(F!=[],
		F=liftall(F);
		if(poldegree(F)==0,
			F=polcoef(F,0););
		res=0;
		for(k=0,poldegree(F,x),
			coef=liftall(polcoef(F,k));
			if((coef!=0)&&(poldegree(coef))!=0,
				print("error:CleanChi: coef not in Q. coef = ",coef," degree coef = ",poldegree(coef));
				);
			res=res+polcoef(coef,0)*x^k;
			);
		print("Chi: polynomial that we get: ",res);
		\\res=polredbest(res);
		\\print("Chi: after polredbest: ",res);
		);	
	return(res);
}

PrepChi(Jinit,CCurve,P0,P1,P2)=
\\WARNING: MODIFY J
{
	my(J,liste);
	J=Jinit;
	J[4]=[Jinit[4][1],Jinit[4][2],[],[]];
	liste=StartingAlpha(Jinit,CCurve,P0,P1,P2);
	J[4][3]=liste[1];
	J[4][4]=liste[2];
	J[14]=liste[3];
	return(J);
}

/* <3 */
InitChiJ2(J,listeRpzOrbiteJ2,listeOrbiteJ2)=
{
	my(listeOrbiteAlpha,Chi,flag,iAlpha,jAlpha,listeCandidat);
	print("beginning: Chi J[2]");
	iAlpha=1;
	jAlpha=2;
	flag=2;
	listeCandidat=vector(2*getd0(J)+1-getg(J),i,i);
	listeCandidat=listeCandidat[3..#listeCandidat];
	while((flag!=1)&&(#listeCandidat!=0),
		print("TentativeChiJ2: iAlpha = ",iAlpha," jAlpha = ",jAlpha);
		[Chi,listeOrbiteAlpha,flag]=TentativeChiJ2(J,listeRpzOrbiteJ2,listeOrbiteJ2,iAlpha,jAlpha);
		print("InitChiJ2: flag = ",flag);
		if(flag==2,
			iAlpha=listeCandidat[1];
			listeCandidat=listeCandidat[2..#listeCandidat];
			);
		if(flag==3,
			jAlpha=listeCandidat[1];
			listeCandidat=listeCandidat[2..#listeCandidat];
			);
		);
	if(#listeCandidat==0,
		print("Warning : error AlphaMap : no combinaisons of (iAlpha,jAlpha) worked");
		);
		
	Chi=CleanChi(Chi);
	if(poldegree(Chi)!=2^(2*getg(J))-1,
		print("WARNING: ChiJ2 error: not the right degree"););
	print("end");
	return([Chi,listeOrbiteAlpha]);
}

InitChiJ2Nicolas(J,PicTorsFrobGendata,Auts)=
\\we just use Nicolas' code to the classic situation but we adapt it to not recompute the data we already computed
{

	my(matAuts,res,res2,listeChiJ2,bestpol,flagChiJ2);
	matAuts=[];
	print("Chi J2 Nicolas: beginning");
	print("Computing the alphas");
	res=PicTorsSpaceFrobEval(J,PicTorsFrobGendata[1],PicTorsFrobGendata[2],2,PicTorsFrobGendata[3],matAuts); \\all the possible alphas
	\\See PicTorsFromBasis fro details
	print("Computing the Chi(s)");
	res2=AllPols(J,res,2,2*getg(J),PicTorsFrobGendata[3]);
	listeChiJ2=vector(#res2,i,res2[i][1]);\\ all the chi
	print("finding the best one");
	flagChiJ2=(#listeChiJ2>0);
	if(flagChiJ2,bestpol=listeChiJ2[1],
		bestpol=0;
		print("error ChiJ2Nicolas: no chi found");
		);
	print("end");
	return([bestpol,flagChiJ2]);
}






\\*********************************************************************** CHI OTC ************************************************************************

InitlisteRpzOrbiteOTCCore(J,listeOrbiteOdd,BasisJ2GAdic,WB0tGAdic,parity)=
\\listeOrbiteOTC sert à la fois aux calculs des alphas et des theta polynomials
{
        my(listeRpzOrbiteOTC,vectornul,Orbite,coef,W);
	if(parity,print("	EVEN TC computation"),
		print("		ODD TC computation"));
        listeRpzOrbiteOTC=vector(#listeOrbiteOdd,i,0);
	vectornul=vector(#BasisJ2GAdic,i,Mod(0,2));
/*
        for(k=1,#listeOrbiteOdd,
                Orbite=listeOrbiteOdd[k];
                coef=Orbite[1];\\print("coef = ",coef);
		if(coef==vectornul,W=WB0tGAdic,
			W=piclc(J,lift(coef),BasisJ2GAdic);
			W=picsum(J,[W,WB0tGAdic]);
			);
		listeRpzOrbiteOTC[k]=W;
		/*
                WCoef=CoefNonZero(coef);print("WCoef = ",WCoef);
                if(WCoef!=[],\\case W!=WB0
                        \\W=SumJ(J,BasisJ2GAdic[WCoef[1]],WB0tGAdic);\\Main difference with TentativeChiJ2
                        W=BasisJ2GAdic[WCoef[1]];
			for(i=2,#WCoef,
                                W=SumJTorsion(J,W,BasisJ2GAdic[WCoef[i]]);
                                );
			W=SumJ(J,W,WB0tGAdic); \\Main difference with TentativeChiJ2
                        listeRpzOrbiteOTC[k]=W;
                        );
		
                if(WCoef==[],\\case W=WB0. WARNING: in opposite with the case J[2], WB0 may not be a rational point. So we need to add it to the polynomial if we want to hope to factor it (over Q) 
                        listeRpzOrbiteOTC[k]=WB0tGAdic;
                        );
		*/
	\\	print("		orbit ",k,"/",#listeOrbiteOdd," done");
         \\       );
	\\A=listeRpzOrbiteOTC;


	parfor(k=1,#listeOrbiteOdd,
                Orbite=listeOrbiteOdd[k];
                coef=Orbite[1];\\print("coef = ",coef);
                if(coef==vectornul,W=WB0tGAdic,
                        W=piclc(J,lift(coef),BasisJ2GAdic);
                        W=picsum(J,[W,WB0tGAdic]);
                        );
		print("         orbit ",k,"/",#listeOrbiteOdd," done");
		W,
			r,listeRpzOrbiteOTC[k]=r);
	
	\\B=listeRpzOrbiteOTC;for(k=1,#A,print("A == B ???",CheckEqual(J,A[k],B[k])));		

        return(listeRpzOrbiteOTC);
}

InitlisteRpzOrbiteOTC(J,listeOrbiteOdd,BasisJ2GAdic,WB0tGAdic,flagTC)=
\\listeOrbiteOTC sert à la fois aux calculs des alphas et des theta polynomials
{
	my(listeRpzOrbiteOTC,listeRpzOrbiteETC);
	if(flagTC,listeRpzOrbiteOTC=InitlisteRpzOrbiteOTCCore(J,listeOrbiteOdd[1],BasisJ2GAdic,WB0tGAdic,0),
		listeRpzOrbiteOTC=InitlisteRpzOrbiteOTCCore(J,listeOrbiteOdd,BasisJ2GAdic,WB0tGAdic,0));
	if(flagTC,
		listeRpzOrbiteETC=InitlisteRpzOrbiteOTCCore(J,listeOrbiteOdd[2],BasisJ2GAdic,WB0tGAdic,1);
		listeRpzOrbiteOTC=[listeRpzOrbiteOTC,listeRpzOrbiteETC];
		);
	return(listeRpzOrbiteOTC);
}







ComputeAlphabyTF(J,RpzTFECoefBasis,iAlpha,jAlpha)=
{
	my(flagTFE,TFECoefBasis,alpha,flagalpha);
	\\[flagTFE,TFECoefBasis]=listeRpzOrbiteOTC[k];
	[flagTFE,TFECoefBasis]=RpzTFECoefBasis;
        if(flagTFE,[alpha,flagalpha]=AlphaMapTF(J,TFECoefBasis,iAlpha,jAlpha),
        	[alpha,flagalpha]=[0,0]);
	return([alpha,flagalpha]);
}

ManageOrbiteAlpha(J,alpha,OrbiteOdd)=
{
	my(Orbitealpha,flag);
	Orbitealpha=OrbiteElt(J,alpha);
        \\listeOrbiteAlpha[k]=Orbitealpha;
        if(default(debug)>2,print("size Orbit alpha = ",#Orbitealpha," = ",#OrbiteOdd));
        \\if(#Orbitealpha!=#listeOrbiteOdd[k],flagglobal=4);
	flag=(#Orbitealpha==#OrbiteOdd);
	if(flag==0,print("ERROR Chi : alphaMap not injective"));
	return([Orbitealpha,flag]);
}

BuildAttemptChi(J,listeOrbiteAlpha,flag)=
{
	my(F);
	F=1;
        if(flag==1,
                for(i=1,#listeOrbiteAlpha,
                        for(j=1,#listeOrbiteAlpha[i],
                                F=F*(x-rpzJtoQq(J,listeOrbiteAlpha[i][j]));
                                );
                        );
                );
	return(F);
}

TentativeChiOTCLoop(J,listeRpzOrbiteOTC,listeOrbiteOdd,iAlpha,jAlpha,flagTF)=
{
        my(listeOrbiteAlpha,flagglobal,W,flagalpha,alpha,Orbitealpha,F,TFECoefBasis,flagOrbite);
        listeOrbiteAlpha=vector(#listeRpzOrbiteOTC,i,0);
        flagglobal=1;

	for(k=1,#listeOrbiteAlpha,

		if(flagglobal==1,
			print("TentativeChiOTCLoop: k = ",k);
			
			if(flagTF==0,
                		W=listeRpzOrbiteOTC[k];
                		[alpha,flagalpha]=AlphaMapTHIBAUT(J,W,iAlpha,jAlpha);
				);
			if(flagTF==1,
				[alpha,flagalpha]=ComputeAlphabyTF(J,listeRpzOrbiteOTC[k],iAlpha,jAlpha);\\WARNING: it's normal: here listeRpzOrbiteOTC=listeRpzTFECoefBasis 
				\\[flagTFE,TFECoefBasis]=listeRpzOrbiteOTC[k];
				\\if(flagTFE,[alpha,flagalpha]=AlphaMapTF(J,TFECoefBasis,iAlpha,jAlpha),
				\\	[alpha,flagalpha]=[0,0]);
				);
			
			if(flagalpha!=1,
                                        flagglobal=flagalpha;);
                	\\alpha=AlphaMapNicolas(J,W,iAlpha,jAlpha);
			if(flagglobal==1,
				[Orbitealpha,flagOrbite]=ManageOrbiteAlpha(J,alpha,listeOrbiteOdd[k]);
				listeOrbiteAlpha[k]=Orbitealpha;
				if(flagOrbite==0,flagglobal=4);
                		\\Orbitealpha=OrbiteElt(J,alpha);
                		\\listeOrbiteAlpha[k]=Orbitealpha;
                		\\print("size Orbit alpha = ",#Orbitealpha," = ",#listeOrbiteOdd[k]);
				\\print("\n");
				\\if(#Orbitealpha!=#listeOrbiteOdd[k],flagglobal=4;);
                		);
			);
		);
	F=BuildAttemptChi(J,listeOrbiteAlpha,flagglobal);
        \\F=1;
	\\if(flagglobal==1,
        \\	for(i=1,#listeOrbiteAlpha,
        \\        	for(j=1,#listeOrbiteAlpha[i],
        \\                	F=F*(x-rpzJtoQq(J,listeOrbiteAlpha[i][j]));
        \\                	);
	\\		);
	\\	);
        \\F=bestappr(F);
        \\IL RESTE A METTRE F DANS LE BON FORMALISME (bestappr, liftall, rendre coef nb et pas pol, polredbest)
        return([F,listeOrbiteAlpha,flagglobal]);
}


TentativeChiOTC(J,listeRpzOrbiteOTC,listeOrbiteOdd,iAlphaI,jAlphaI,maxindex,largestindex,flagTF)=
{
	my(iAlpha,jAlpha,flag,listeCandidat,Chi,listeOrbiteAlpha);
	iAlpha=iAlphaI;
	jAlpha=jAlphaI;
	flag=2;
        listeCandidat=vector(largestindex,i,i);
        listeCandidat=listeCandidat[(maxindex+1)..#listeCandidat];
        while((flag!=1)&&(#listeCandidat!=0),
                print("TentativeChiOTC: iAlpha = ",iAlpha," jAlpha = ",jAlpha);
                [Chi,listeOrbiteAlpha,flag]=TentativeChiOTCLoop(J,listeRpzOrbiteOTC,listeOrbiteOdd,iAlpha,jAlpha,flagTF);print("InitChiOTC: flag = ",flag);
                if(flag==4,
			iAlpha=listeCandidat[1];
                        listeCandidat=listeCandidat[2..#listeCandidat];
			print("Map Alpha not injective");
			);
		if(flag==2,
                        iAlpha=listeCandidat[1];
                        listeCandidat=listeCandidat[2..#listeCandidat];
			print("Map Alpha not defined");
                        );
                if(flag==3,
                        jAlpha=listeCandidat[1];
                        listeCandidat=listeCandidat[2..#listeCandidat];
			print("Map Alpha not defined");
                        );
                );
        if(#listeCandidat==0,
                print("Warning : error AlphaMap : no combinaisons of (iAlpha,jAlpha) worked");
		flag=4;
                );
	return([flag,Chi,listeOrbiteAlpha,iAlpha,jAlpha]);
}

PolishChiTC(g,ChiInit,parity)=
{
	my(Chi,CardOTC,flagdeg,flagroot,flag);
	Chi=CleanChi(ChiInit);
	if(parity==1,CardOTC=(2^(g-1))*(2^g-1),
		CardOTC=(2^(g-1))*(2^g+1));
        flagdeg=(poldegree(Chi)!=CardOTC);
                \\flagroot=(gcd(Chi,Chi')!=1);
        flagroot=(poldegree(gcd(Chi,Chi'))!=0);
	if(flagdeg||flagroot,print("flagdeg = ",flagdeg," flagroot = ",flagroot," poldegree = ",poldegree(Chi)," CardOTC = ",CardOTC," parity = ",parity," gcd( Chi, Chi')"," = ",gcd(Chi,Chi')," poldegree(gcd(Chi,Chi'): ",poldegree(gcd(Chi,Chi')) ));
	flag=(flagdeg==0)&&(flagroot==0);
	return([Chi,flagdeg,flagroot,flag]);
}

PolishChi(g,ChiInit)={return(PolishChiTC(g,ChiInit,1));}

/* <3 */
InitChiOTC(J,listeRpzOrbiteOTC,listeOrbiteOdd,flagTF)=
\\if flagTF=1, we build ChiOTC using the Theta Fractions
\\		BECAREFULL: In that case, listeRpzOrbiteOTC=listeRpzTFECoefBasis
\\if flagTF=0, we use the classical method
{
        my(iAlpha,jAlpha,maxindex,largestindex,flag,listeCandidat,listeOrbiteAlpha,Chi,flagdeg,flagroot);
	if(default(debug)>1,print("beginning: chi OTC"));
	iAlpha=1;
        jAlpha=2;
	maxindex=max(iAlpha,jAlpha);
	if(flagTF==0,largestindex=2*getd0(J)+1-getg(J),
		largestindex=getg(J));
	flag=0;
	while((flag!=1)&&(maxindex<largestindex),
		[flag,Chi,listeOrbiteAlpha,iAlpha,jAlpha]=TentativeChiOTC(J,listeRpzOrbiteOTC,listeOrbiteOdd,iAlpha,jAlpha,maxindex,largestindex,flagTF);
		maxindex=max(iAlpha,jAlpha);

		[Chi,flagdeg,flagroot,flag]=PolishChi(g,Chi);	
        	\\Chi=CleanChi(Chi);
		\\flagdeg=(poldegree(Chi)!=(2^(getg(J)-1))*(2^getg(J)-1));
		\\flagroot=(gcd(Chi,Chi')!=1);
		\\flagroot=(poldegree(gcd(Chi,Chi'))!=0);
		if(flagdeg||flagroot,
		\\if(poldegree(Chi)!=(2^(getg(J)-1))*(2^getg(J)-1),
			\\print("poldegree = ",poldegree(Chi)," gcd( Chi, Chi')"," = ",gcd(Chi,Chi')); 
			flag=0;
			if(iAlpha<jAlpha,
				iAlpha=iAlpha+1;
				if(iAlpha==jAlpha,
					jAlpha=jAlpha+1;
					);
				);
			if(iAlpha>jAlpha,
				jAlpha=jAlpha+1;
				if(iAlpha==jAlpha,
					iAlpha=iAlpha+1;
					);
				);
			maxindex=maxindex+1;			
			print("WARNING: ChiOTC error: ChiOTC is not of the right degree (Interpolation failed)");
			);
		);
	if(flag,PrintOrbitDecompositionChi(Chi));
	print("end");
        return([Chi,listeOrbiteAlpha,flag]);
}













InitClassicChiOTC(J,listeRpzOrbiteOTC,listeOrbiteOdd)=
{
	print("\n");
	print("Chi OTC Classic");
	return( InitChiOTC(J,listeRpzOrbiteOTC,listeOrbiteOdd,0) );
}












CleanAllPols(res2)=
{
	my(listeChi,bestpol,flagglobal,listeRpzAlpha);
	listeChi=vector(#res2,i,res2[i][1]);\\ all the chi
        print("finding the best one");
	flagglobal=(#listeChi>0);
        if(flagglobal,bestpol=listeChi[1],
                bestpol=0;
                print("error InitChiOTCNicolas: no chi found");
                );
	listeRpzAlpha=0;
	if(flagglobal,
		listeRpzAlpha=res2[1][2]);
	return([bestpol,flagglobal,listeRpzAlpha]);
}

RecoverlisteOrbiteAlpha(J,listeRpzAlpha,listeOrbite,flagglobalInit)=
{
	my(flagglobal,listeOrbiteAlpha,Orbitealpha,flagOrbite,alpha,index);
	flagglobal=flagglobalInit;
        listeOrbiteAlpha=vector(#listeOrbite,i,0);
	index=1;
        for(k=1,#listeOrbite,
		alpha=listeRpzAlpha[index];
                \\[Orbitealpha,flagOrbite]=ManageOrbiteAlpha(J,listeAlpha[k],listeOrbiteOdd[k]);
                [Orbitealpha,flagOrbite]=ManageOrbiteAlpha(J,alpha,listeOrbite[k]);
		listeOrbiteAlpha[k]=Orbitealpha;
		index=index+#Orbitealpha;
                if(flagOrbite==0,flagglobal=4);
                );
	return([listeOrbiteAlpha,flagglobal]);
}


CheckAlpha(listeAlpha)=
{
	my(flag);
	flag=1;
	for(k=1,#listeAlpha,
		if(listeAlpha[k]==0,
			flag=0;));
	return(flag);
}

AlphaGood(J,listeRpzOrbite,CCurve,P0,P1,P2,Frob)=
{
	my(listeAlpha);
	\\J=JInit;
	\\J=PrepChi(J,CCurve,P0,P1,P2,Frob);
	\\listeAlpha=parapply(W->PicEval(J,W),listeRpzOrbite);
	listeAlpha=vector(#listeRpzOrbite,i,0);
	
	for(i=1,#listeRpzOrbite,
		print("new alpha: ",i,"/",#listeRpzOrbite);
                listeAlpha[i]=PicEval(J,listeRpzOrbite[i]);
		);

	/*
	parfor(i=1,#listeRpzOrbite,
		print("new alpha: ",i,"/",#listeRpzOrbite);
		PicEval(J,listeRpzOrbite[i]),
			r,listeAlpha[i]=r);
	*/

	if((CheckAlpha(listeAlpha)==0),print("\n Error Alpha \n"));
	/*
	if((CheckAlpha(listeAlpha)==0)&&(P2!=0)&&(P1!=0),
		print("trying another WE1/WE2");
                J=PrepChi(J,CCurve,P0,P1,0,Frob);
		listeAlpha=parapply(W->PicEval(J,W),listeRpzOrbite);
		);
	if((CheckAlpha(listeAlpha)==0)&&(P2!=0),
		print("trying another WE1/WE2");
                J=PrepChi(J,CCurve,P0,P2,0,Frob);
		listeAlpha=parapply(W->PicEval(J,W),listeRpzOrbite);
                );
	if((CheckAlpha(listeAlpha)==0)&&(P1!=0),
		print("trying another WE1/WE2");
                J=PrepChi(J,CCurve,P0,0,0,Frob);
		listeAlpha=parapply(W->PicEval(J,W),listeRpzOrbite);
                );
	*/
	return([J,CheckAlpha(listeAlpha),listeAlpha]);
}

InitChiNicolasCore(JInit,listeRpzOrbite,listeOrbite,CCurve,P0,P1,P2,Frob)=
\\we just use Nicolas' code to the classic situation but we adapt it to not recompute the data we already computed
\\listeOrbite is liste Orbite TC
{
	my(J,res,res2,bestpol,listSizeOrbit,listeAlphaVec,flagglobal,listeOrbiteAlpha,listeRpzAlpha);
        print("Chi OTC Nicolas: beginning");
        print("Computing the alphas");
        \\res=PicTorsSpaceFrobEval(J,PicTorsFrobGendata[1],PicTorsFrobGendata[2],2,PicTorsFrobGendata[3],matAuts); \\all the possible alphas
       
	\\nmodF=#listeRpzOrbiteOTC;
	listSizeOrbit=vector(#listeOrbite,i,#listeOrbite[i]);
	\\ImodF=listSizeOrbit;
	\\res=AlphaNicolas(nmodF,[J],listeRpzOrbiteOTC,ImodF);
	\\ATTENTION: there still is the problem with ImodF2:=res[2]
	J=JInit;
	[J,flagglobal,listeAlphaVec]=AlphaGood(J,listeRpzOrbite,CCurve,P0,P1,P2,Frob);
	if(flagglobal==0,
		print("\n","ERROR alpha","\n"));
	res=[listeAlphaVec,Vecsmall(listSizeOrbit)];

	print("Computing the Chi(s)");
        res2=AllPols(J,res,2,0,0);

	[bestpol,flagglobal,listeRpzAlpha]=CleanAllPols(res2);
        print("Chi Nicolas: ",bestpol);
        print("GOOD CHI");
        PrintOrbitDecompositionChi(bestpol);

	[listeOrbiteAlpha,flagglobal]=RecoverlisteOrbiteAlpha(J,listeRpzAlpha,listeOrbite,flagglobal);

        print("end");

	\\return([ChiOTC,listeOrbiteAlpha,flagChiOTC]) 
        return([bestpol,listeOrbiteAlpha,flagglobal]);
}

InitChiOTCNicolas(J,listeRpzOrbiteOTC,listeOrbiteOdd,CCurve,P0,P1,P2,Frob)=
{
	return(InitChiNicolasCore(J,listeRpzOrbiteOTC,listeOrbiteOdd,CCurve,P0,P1,P2,Frob));
}















\\*********************************************** LAMBDA *******************************************************************

export(InitLambdaBetaGAdic)
export(CheckBetaEffective)
export(RedJ)
export(RedPointJ)
export(DivSubDimKFF)
export(KtoWFF)
export(reduceFFMatrix)
export(liftFFMatrix)
export(getWE)
export(getWn)
export(getd0)
export(getg)
export(mulPointWise)
export(WtoKFF)
export(reduceFF)
export(liftFF)
export(MulFF)
export(RpzLambdaParforAux)
export(Choicewarf)
export(RpzLambdaParfor)

Choicewarf(WArfGAdic,dimArf,flagAlpha,parity)=
\\warf is for OTC or vanishing ETC
\\warf2 is for the non-vanishing case
\\we adjust in fonction of flagAlpha, i.e. normal Alpha or Alpha with one more condition on P0
{
	my(warf,warf2,l,dimarf,dimarf2);
	if(type(WArfGAdic)=="t_MAT",WArfGAdic=[WArfGAdic,0,0,0]);
        if(parity,
        	if(flagAlpha,l=1,l=3);
                warf=WArfGAdic[l];
		dimarf=dimArf[l];
                warf2=0;
		dimarf2=-1;
                );
	if(parity==0,
        	if(flagAlpha,l=2,l=1);
                warf=WArfGAdic[l];\\vanishing case
		dimarf=dimArf[l];
                if(flagAlpha,l=3,l=4);
                warf2=WArfGAdic[l];\\non vanishing case
		dimarf2=dimArf[l];
                );
	dimarf=[dimarf,dimarf2];
	return([warf,warf2,dimarf]);
}

RpzLambdaParforAux(J,WTC,warf,warf2,dimarf,parity,flagETCnonVanish,Gx,flagAlpha,W2plusP0)=
{
	my(dimBeta,LambdaGAdic,flagWE1,flagLambda);
	[dimBeta,LambdaGAdic,flagWE1]=InitLambdaBetaGAdic(J,WTC,warf,dimarf[1],parity,0,Gx,flagAlpha,W2plusP0);
        if(default(debug)>2,print("dim Beta = ",dimBeta));
        if((parity==0)&&flagETCnonVanish&&(dimBeta!=1),
        	print("\n","Non Vanishing Beta: we try to impose P0 as a zero","\n");
                [dimBeta,LambdaGAdic,flagWE1]=InitLambdaBetaGAdic(J,WTC,warf2,dimarf[2],parity,1,Gx,flagAlpha,W2plusP0);\\vanishing=1
                );
	flagLambda=(dimBeta==1)&&flagWE1;
        if(dimBeta!=1,
        	flagLambda=0;
                LambdaGAdic=0;
                );
	return([flagLambda,LambdaGAdic]);
}

/*
RpzLambdaParfor(J,listeRpzOrbiteOTC,WArfGAdic,dimArf,parity,flagETCnonVanish,W2plusP0,listeRpzGx,i)=
{
	my(flagAlpha,warf,warf2,dimarf,Gx,Lambda);
        flagAlpha=listeRpzGx[i][1];
        print("flagAlpha = ",flagAlpha);
        [warf,warf2,dimarf]=Choicewarf(WArfGAdic,dimArf,flagAlpha,parity);
        Gx=listeRpzGx[i][2];
        Lambda=RpzLambdaParforAux(J,listeRpzOrbiteOTC[i],warf,warf2,dimarf,parity,flagETCnonVanish,Gx,flagAlpha,W2plusP0);
	print("Lamba: GAdic orbit (",i,"/",#listeRpzOrbiteOTC,") done");
	return(Lambda);
}
*/

InitlisteRpzLambdaGAdicCore(J,listeRpzOrbiteOTC,WArfGAdic,dimArf,parity,flagETCnonVanish,listeRpzGx,W2plusP0)=
{
	my(listeRpzLambdaGAdic,flagtotal,flagAlpha,Gx,warf,warf2,dimarf);
        if(parity==0,print("	EVEN TC computation"),print("	ODD TC computation"));
	listeRpzLambdaGAdic=vector(#listeRpzOrbiteOTC,i,0);

/*
        for(i=1,#listeRpzOrbiteOTC,
                WTC=listeRpzOrbiteOTC[i];
                print("Lamba: new GAdic orbit (",i,"/",#listeRpzOrbiteOTC,")");
		if(parity,warf=WArfGAdic,
			if(flagETCnonVanish==0,warf=WArfGAdic,
				warf=WArfGAdic[1]);
			);
		\\print("warf: size: ",type(warf));
                [dimBeta,LambdaGAdic,flagWE1]=InitLambdaBetaGAdic(J,WTC,warf,parity,0);
                print("dim Beta = ",dimBeta);
		if((parity==0)&&flagETCnonVanish&&(dimBeta!=1),
			print("\n","Non Vanishing Beta: we try to impose P0 as a zero","\n");
			warf=WArfGAdic[2];
			[dimBeta,LambdaGAdic,flagWE1]=InitLambdaBetaGAdic(J,WTC,warf,parity,1);\\vanishing=1
			);
		flagLambda=(dimBeta==1)&&flagWE1;
                if(dimBeta!=1,
                        flagLambda=0;
                        LambdaGAdic=0;
                        );
                listeRpzLambdaGAdic[i]=[flagLambda,LambdaGAdic];
                );
	A=listeRpzLambdaGAdic;
*/






	if(default(debug)>1,print("why 'parfor' is not working here and I have to use a normal loop 'for' ?"));
	for(i=1,#listeRpzOrbiteOTC,
		print("new representative (",i,"/",#listeRpzOrbiteOTC,")");
		flagAlpha=listeRpzGx[i][1];
		if(default(debug)>1,print("flagAlpha = ",flagAlpha));
		[warf,warf2,dimarf]=Choicewarf(WArfGAdic,dimArf,flagAlpha,parity);
		Gx=listeRpzGx[i][2];
		listeRpzLambdaGAdic[i]=RpzLambdaParforAux(J,listeRpzOrbiteOTC[i],warf,warf2,dimarf,parity,flagETCnonVanish,Gx,flagAlpha,W2plusP0);
		);
/*
	\\A=listeRpzLambdaGAdic;
	parfor(i=1,#listeRpzOrbiteOTC,
                print("Lamba: new GAdic orbit (",i,"/",#listeRpzOrbiteOTC,")");
                flagAlpha=listeRpzGx[i][1];
                [warf,warf2,dimarf]=Choicewarf(WArfGAdic,dimArf,flagAlpha,parity);
                Gx=listeRpzGx[i][2];
                RpzLambdaParforAux(J,listeRpzOrbiteOTC[i],warf,warf2,dimarf,parity,flagETCnonVanish,Gx,flagAlpha,W2plusP0),
			r,listeRpzLambdaGAdic[i]=r);
	\\B=listeRpzLambdaGAdic;print("A == B ??",A==B);
*/

/*
	\\BECAREFUL: use the function RpzLambdaParfor that has been made a comment. Please uncomment it to use the following.
	\\for(i=1,#listeRpzOrbiteOTC,listeRpzLambdaGAdic[i]=RpzLambdaParfor(J,listeRpzOrbiteOTC,WArfGAdic,dimArf,parity,flagETCnonVanish,W2plusP0,listeRpzGx,i));
	parfor(i=1,#listeRpzOrbiteOTC,
		RpzLambdaParfor(J,listeRpzOrbiteOTC,WArfGAdic,dimArf,parity,flagETCnonVanish,W2plusP0,listeRpzGx,i),
			r,listeRpzLambdaGAdic[i]=r);
*/
	
	flagtotal=1;
	for(i=1,#listeRpzOrbiteOTC,
		if(listeRpzLambdaGAdic[i][1]==0,flagtotal=0));
        return([listeRpzLambdaGAdic,flagtotal]);
}

InitlisteRpzLambdaGAdic(J,listeRpzOrbiteOTC,WArfGAdic,dimArf,flagTC,flagETCeffective,parity,flagETCnonVanish,listeRpzGx,W2plusP0)=
{
	my(listeRpzLambdaGAdic,flag);

	if(flagTC==0,
		[listeRpzLambdaGAdic,flag]=InitlisteRpzLambdaGAdicCore(J,listeRpzOrbiteOTC,WArfGAdic,dimArf,1,0,listeRpzGx,W2plusP0));

	if(flagTC,
		if(parity==1,i=1);
		if(parity==0,i=2);
		/*
		if(flagETCeffective==0,warf=WArfGAdic,
			warf=WArfGAdic[i];
			if(flagETCnonVanish&&(parity==0),
				warf=[warf,WArfGAdic[3]]);
			);
		*/
		[listeRpzLambdaGAdic,flag]=InitlisteRpzLambdaGAdicCore(J,listeRpzOrbiteOTC[i],WArfGAdic,dimArf,parity,flagETCnonVanish,listeRpzGx[i],W2plusP0);
		);
	return([listeRpzLambdaGAdic,flag]);
}



InversesThetaInvMakdisi(J,sThetaInv)=
{
        my(wrongindex,sTheta);
        wrongindex=[];
        sTheta=vector(#sThetaInv,k,0)~;
        for(k=1,#sTheta,
                \\print("k= ",k);
                \\we test if sThetaInv[k] is inversible. So we look if its reduction is 0. Since p is of good reduction, we check p/lift(a)
                if(Mod(sThetaInv[k],getpMak(J))!=0,sTheta[k]=padicinvElt(J,sThetaInv[k]), \\sTheta[k]=divpadicElt(J,1,sThetaInv[k]),
                        wrongindex=concat(wrongindex,[k]);
                        sTheta[k]=+oo;
                        );
                );
        return([sTheta,wrongindex]);
}

InitTheta(J,ThetaInv)=
{
        my(wrongindex,sTheta);
        [sTheta,wrongindex]=InversesThetaInvMakdisi(J,ThetaInv);
        return([sTheta,wrongindex]);
}

InitTFEGAdic(J,wrongindex,lambda,sTheta)=
{
        my(TFEGAdic,c);
        c=vector(#sTheta,l,0)~;
        for(l=1,#c,
        	if(IsInList(wrongindex,l),c[l]=+oo,
                	\\print("sTheta == 1 ? ",sTheta[l]==sTheta[1]);
                        c[l]=mulpadicElt(J,sTheta[l],mulpadicElt(J,lambda[l],lambda[l])); \\f*lambda^2
                                \\print("lambda nul ? ",lambda[l]==0);
                                \\print("c == 1 ? ",c[l]==c[1]);
                        );
                );
	\\TFEGAdic=matrix(#c,1,i,j,c[i]);
        TFEGAdic=c;
	return(TFEGAdic);
}

InitTFECoefBasis(J,TFEGAdic,WK0,wrongindex)=
{
        my(TFECoefBasis,flagwi);
        [flagwi,TFECoefBasis]=ColtoPolCoefClean(J,TFEGAdic,WK0,wrongindex);
        return([flagwi,TFECoefBasis]);
}

InitlisteRpzTFECoefBasisCore(J,listeRpzLambdaGAdic,listeRpzGx,listeNum,WK0)=
{
        my(listeRpzTFECoefBasis,flagLambda,LambdaGAdic,flagTFE,TFECoefBasis,ThetaInv,sTheta,wrongindex,TFEGAdic,flagwi,flagAlpha,Gx);
	if(default(debug)>1,print("ColtoPolCoefGAdicH : we don't assess inversibility of the pivot, only non-nullity"));
        listeRpzTFECoefBasis=vector(#listeRpzLambdaGAdic,i,0);
        for(i=1,#listeRpzLambdaGAdic,
		if(default(debug)>2,print("i= ",i));
                [flagLambda,LambdaGAdic]=listeRpzLambdaGAdic[i];
                flagTFE=1;
                if(flagLambda==0,
                        print("dim Beta != 1. orbit #",i,". We can't compute ChiOTC by the TF. Indeed f_Beta Galois faithfully represents Beta but if dim(Beta)!=1 we can't ensure that Beta is defined over Qbar");
                        flagTFE=0;
                        TFECoefBasis=[];
                        );
                if(flagLambda==1,
			flagAlpha=listeRpzGx[i][1];
			Gx=listeRpzGx[i][2];
			if(flagAlpha==0,print("flagAlpha = 0, it may bug here ??"));
                        ThetaInv=InitsThetaDivSubInv(J,Gx,listeNum[i]);
                        [sTheta,wrongindex]=InitTheta(J,ThetaInv);
                        \\OrbiteTFEGAdic=InitOrbiteTFEGAdic(J,wrongindex,OrbiteLambdaGAdic,lOrbsTheta);\\print("c == 1 ? ",OrbiteTFEGAdic[1][1]==OrbiteTFEGAdic[1][2]);
                        TFEGAdic=InitTFEGAdic(J,wrongindex,LambdaGAdic,sTheta);
                        \\[flagwi,OrbiteTFECoefBasis]=InitOrbiteTFECoefBasis(J,OrbiteTFEGAdic,WK0,wrongindex);\\print("CoefBasis = ",OrbiteTFECoefBasis[1]);
                        [flagwi,TFECoefBasis]=InitTFECoefBasis(J,TFEGAdic,WK0,wrongindex);
                        if(flagwi==0,
                                print("we couldn't identify the TFE because too much poles on the Makdisi evaluation Points Z. orbit #",i);
                                flagTFE=0;
                                TFECoefBasis=[];
                                );
                        );
                listeRpzTFECoefBasis[i]=[flagTFE,TFECoefBasis];
                );
        return(listeRpzTFECoefBasis);
}

InitlisteRpzTFECoefBasis(J,listeRpzLambdaGAdic,listeRpzGx,listeNum,VTFE,flagTC,flagETCeffective)=
\\VTFE=WK0 in case of OTC
{
	my(listeRpzTFECoefBasis,lETC,WK0);
	if(flagETCeffective==0,
		if(flagTC,listeRpzTFECoefBasis=InitlisteRpzTFECoefBasisCore(J,listeRpzLambdaGAdic,listeRpzGx[1],listeNum[1],VTFE),
			listeRpzTFECoefBasis=InitlisteRpzTFECoefBasisCore(J,listeRpzLambdaGAdic,listeRpzGx,listeNum,VTFE));
		);

	if(flagETCeffective,
		if(flagTC==0,print("ERROR: flagTC should be 1"));
		listeRpzTFECoefBasis=InitlisteRpzTFECoefBasisCore(J,listeRpzLambdaGAdic[1],listeRpzGx[1],listeNum[1],VTFE[1]);
		lETC=InitlisteRpzTFECoefBasisCore(J,listeRpzLambdaGAdic[2],listeRpzGx[2],listeNum[2],VTFE[2]);
		listeRpzTFECoefBasis=[listeRpzTFECoefBasis,lETC];
		);
	return(listeRpzTFECoefBasis);
}


AlphaMapTF(J,TFECoefBasis,iAlpha,jAlpha)=
{
	my(num,denom,flagalpha,alpha);
	num=TFECoefBasis[iAlpha];
	denom=TFECoefBasis[jAlpha];
	flagalpha=(Mod(denom,getpMak(J))!=0);
	if(flagalpha,alpha=divpadicElt(J,num,denom),
		alpha=0);
	return([alpha,flagalpha]);
}

AlphaMapTFCLCore(J,TFECoefBasis,lf)=
{
	my(a);
	a=0;for(k=1,#lf,a=a+mulpadicElt(J,lf[k],TFECoefBasis[k]));
	return(a);
}

AlphaMapTFCL(J,TFECoefBasis,lf1,lf2)=
{
	my(num,denom,alpha);
	num=AlphaMapTFCLCore(J,TFECoefBasis,lf1);
	denom=AlphaMapTFCLCore(J,TFECoefBasis,lf2);
	alpha=divpadicElt(J,num,denom);
	return(alpha);
}

InitChiOTCbyTF(J,listeRpzTFECoefBasis,listeOrbiteOdd)=
{
	print("\n");
	print("Chi (O)TC by Theta Fractions");
	return(InitChiOTC(J,listeRpzTFECoefBasis,listeOrbiteOdd,1));
}












randomQ(bound)=
{
        \\return((random(bound)-random(bound))/(1+random(bound/2)));
	return(random(bound)-random(bound));
}

CheckNoCoefZero(v,p,emax)=
\\we want all the coefficients are inversible
{
        my(flag);
        flag=1;
        for(k=1,#v,
                if(Mod(v[k]*(1+O(p^emax)),p)==0,flag=0));
        return(flag);
}

CheckColinearlf(lf,lfbis)=
{
        my(flag,index,quotient);
        if(lfbis==0,flag=0,flag=1);
	if(flag,
        	index=FirstNonZero(lf);
        	for(k=1,index-1,if(lfbis[k]!=0,flag=0));
		);
        if(flag,
                quotient=lf[index]/lf[index];
                for(k=index+1,#lf,
                if(flag,
                        if(lf[k]==0,flag=(lfbis[k]==0),
                                flag=((lfbis[k]/lf[k])==quotient)));
                        );
                );
        return(flag);
}

RandomlfBasic(g)=
{
	my(lf);
	lf=vector(g,i,0);
	lf[random(g)+1]=1;
	return(lf);
}

RandomlfSimple(g)=
{
	return(vector(g,i,random(2)-random(2)));
}

RandomlfGeneral(g,bound)=
{
	return(vector(g,i,randomQ(bound)));
}

NonTriviallfCore(J,Malpha,lf1,char)=
\\if lf1!=0, we avoid it
{
        my(countermax,counter,nblign,p,v,flaglf,lf);
        countermax=50;
        counter=0;
        flaglf=0;
	nblign=#Malpha[,1];
	p=getpMak(J);
        while((counter<countermax)&&(flaglf==0),
		if(char=="Basic",
			lf=RandomlfBasic(nblign);
			);
		if(char=="Simple",
			lf=RandomlfSimple(nblign);
			);
		if(char=="General",
                	\\lf=vector(getg(J),i,randomQ(p));
			lf=RandomlfGeneral(nblign,p);
			);
                counter=counter+1;\\print("counter lf = ",counter);
		v=vector(#Malpha,j,AlphaMapTFCLCore(J,Malpha[,j],lf));
                flaglf=CheckNoCoefZero(v,p,geteJ(J))&&(CheckColinearlf(lf,lf1)==0);
                );
	\\print("counter lf = ",counter);
	if(counter==countermax,flaglf=0);
        return([flaglf,lf]);
}

NonTriviallfLoop(J,Malpha,char)=
{
	my(flaglf,lf1,flag,lf2);
	flag=1;
	[flaglf,lf1]=NonTriviallfCore(J,Malpha,0,char);
        if(flaglf==0,flag=0,
        	[flaglf,lf2]=NonTriviallfCore(J,Malpha,lf1,char);
                if(flaglf==0,flag=0);
                );
	return([flag,lf1,lf2]);
}

NonTriviallf(J,Malpha,char)=
{
	my(flaglf,lf1,flag,lf2,countermax);
        print("Linear Forms for projection Nu: ",char);
	flag=0;
	countermax=10;
	for(k=1,countermax,
		if(flag==0,
			[flag,lf1,lf2]=NonTriviallfLoop(J,Malpha,char);
			);
		);
	return([flag,lf1,lf2]);
}


InitChiOTCbyTFwithCLCore(J,Malpha,listeOrbiteOdd,char,parity)=
{
        my(flag,nblign,dim,flaglf,lf1,lf2,listeOrbiteAlpha,alpha,Orbitealpha,Chi,flagdeg,flagroot);
	flag=1;
        nblign=#Malpha[,1];
        dim=#Malpha;
\\Linear Forms
        if(flag,
		[flag,lf1,lf2]=NonTriviallf(J,Malpha,char);
		/*	
		[flag,lf1,lf2]=NonTriviallf(J,Malpha,"Basic");
		if(flag==0,
                	[flag,lf1,lf2]=NonTriviallf(J,Malpha,"Simple");
                	);
		if(flag==0,
                	[flag,lf1,lf2]=NonTriviallf(J,Malpha,"General");
                	);
		*/
		);
	\\print("got lf");
\\Alphas
	if(flag,
		listeOrbiteAlpha=vector(dim,i,0);
		\\print("flag = ",flag);
		for(k=1,dim,
			if(flag,
				\\alpha=(lf1*Malpha[,k])/(lf2*Malpha[,k]);
				alpha=AlphaMapTFCL(J,Malpha[,k],lf1,lf2);
				[Orbitealpha,flag]=ManageOrbiteAlpha(J,alpha,listeOrbiteOdd[k]);
				if(flag,listeOrbiteAlpha[k]=Orbitealpha);
				);
			);
		);
	\\print("got listeOrbiteAlpha");
\\Chi
	\\print("flag = ",flag);
	if(flag,
		Chi=BuildAttemptChi(J,listeOrbiteAlpha,flag);\\print("Chi = ",Chi);
		[Chi,flagdeg,flagroot,flag]=PolishChiTC(getg(J),Chi,parity);\\print("Chi = ",Chi);
		\\flag=(flagdeg==0)&&(flagroot==0);
		);
	if(flag==0,
		Chi=0;
		listeOrbiteAlpha=vector(dim,i,0);
		);
        return([Chi,listeOrbiteAlpha,flag]);
}

ChiOTCbyTFwithCLAttempt(J,Malpha,listeOrbiteOdd,char,parity)=
{
	my(countermax,counter,flag,Chi,listeOrbiteAlpha);
	countermax=5;
        counter=0;
        flag=0;
        while((flag==0)&&(counter<countermax),
        	counter=counter+1;
                [Chi,listeOrbiteAlpha,flag]=InitChiOTCbyTFwithCLCore(J,Malpha,listeOrbiteOdd,char,parity);
                );
	if(default(debug)>2,print("lf attempt: ",counter));
	return([Chi,listeOrbiteAlpha,flag]);
}

InitChiOTCbyTFwithCL(J,listeRpzTFECoefBasis,listeOrbiteOdd,parity)=
{
	my(nblign,dim,Malpha,flag,countermax,counter,Chi,lOrbiteAlpha);
	if(parity==0,print("	EVEN TC computation"),print("	ODD TC computation"));
        nblign=getg(J);
	for(k=1,#listeRpzTFECoefBasis,
		if(type(listeRpzTFECoefBasis[k])=="t_VEC",
			nblign=#listeRpzTFECoefBasis[k][2]));
	flag=1;
        dim=#listeRpzTFECoefBasis;
\\Malpha
        Malpha=matrix(nblign,dim,i,j,0);
        for(k=1,dim,
                if(listeRpzTFECoefBasis[k][1],Malpha[,k]=listeRpzTFECoefBasis[k][2]~,
                        print("ERROR: TFECoefBasis wrong");
                        flag=0;
                        );
                );
	if(flag,
		[Chi,lOrbiteAlpha,flag]=ChiOTCbyTFwithCLAttempt(J,Malpha,listeOrbiteOdd,"Basic",parity);
		if(flag==0,
			[Chi,lOrbiteAlpha,flag]=ChiOTCbyTFwithCLAttempt(J,Malpha,listeOrbiteOdd,"Simple",parity);
			);
		if(flag==0,
                        [Chi,lOrbiteAlpha,flag]=ChiOTCbyTFwithCLAttempt(J,Malpha,listeOrbiteOdd,"General",parity);
                        );
		/*
		countermax=10;
		counter=0;
		flag=0;
		while((flag==0)&&(counter<countermax),
			counter=counter+1;
			[Chi,listeOrbiteAlpha,flag]=InitChiOTCbyTFwithCLCore(J,Malpha,listeOrbiteOdd);
			);
		*/
		);
	if(flag,print("GOOD Chi"),
                Chi=0;
                lOrbiteAlpha=vector(dim,i,0);
                );
	if(flag,PrintOrbitDecompositionChi(Chi));
	print("end");
	return([Chi,lOrbiteAlpha,flag]);
}





















\\*********************************************************************************************************************************************
\\*********************************************************************************************************************************************










\\TEST
TestWB0tGAdicK0(J,WB0tGAdic,WK0tGAdic)=
{
	return(CheckEqual(J,MulInt(J,WB0tGAdic,2),WK0tGAdic));
}

TestWB0tGAdicRed(J,WB0tGAdicHigh,WB0tGAdicLow,eLow)=
{
	my(Jold);
	Jold=RedJ(J,eLow);
	return(CheckEqual(Jold,WB0tGAdicHigh,WB0tGAdicLow));
}
\\






PrepChiJ2(J,Frob,SBasisJ2GAdic)=
{
	my(listeOrbiteJ2,listeRpzOrbiteJ2);
	print("Orbite J[2]: rpz vector form");
        listeOrbiteJ2=InitlisteOrbiteJ2(getg(J),Frob);
        print("end");
        print("Orbite J[2]: Makdisi form (long and useful only for ChiJ[2])");
        listeRpzOrbiteJ2=InitlisteRpzOrbiteJ2(J,listeOrbiteJ2,SBasisJ2GAdic); \\(long and useful only for ChiJ[2])
        print("end");
	return([listeRpzOrbiteJ2,listeOrbiteJ2]);
}








\\******************* MAIN ******************************************************************
\\*******************************************************************************************


\\PrepInitOTCGAdic(J,JFF,WB0tFF,SBasisJ2FF,P0,CCurve,CCurvep,WK0tGAdic,CardJq,CardJp)=
	\\my(Jfinal,g,fpol,Frob,A,M,S,v0,testArf,ArfB0,WeilMatrix,WArf,Arfdata,B0data,flagArf,listeOrbiteOdd,listOTCNt,listeRpzOrbiteOTC,WB0tGAdic,SBasisJ2GAdic);
	
DataOTCLocal(JFF,CCurvep,P0,SBasisJ2FF,WB0tFF,flagTC)=
{
	my(g,fpol,WeilMatrix,WArf,Arfdata,Frob,A,M,S,ArfB0,v0,testArf,B0data,flagArf,listeOrbiteOdd,listOTCNt);
	if(default(debug)>1,print("\n","beginning: OTC")); 
	g=getg(JFF);
	fpol=Mod(JFF[1],getpMak(JFF));
        \\CCurvep=crvinit(fpol,'t,'a);
	WeilMatrix=InitSymplecticMatrix(g);
	if(default(debug)>1,print("beginning InitWArf"));
	WArf=InitWArf(JFF,CCurvep,P0);
	if(default(debug)>1,print("end"));
	Arfdata=[WArf,CCurvep,P0];
	if(default(debug)>1,print("beginning Frob Matrix"));
	\\Frob=InitFrob(JFF,g,SBasisJ2FF,Weildata,WeilMatrix,CCurvep);
	Frob=InitFrob(JFF,g,SBasisJ2FF,WeilMatrix);
	if(default(debug)>1,print("end"));
	A=InitMatA(Frob);\\A = Id-Frob
	M=matimage(A);
	S=SuppImA(M,A);

	print("finding the dictionnary TC <-> QTF ");
	\\ArfB0=QuadCanonique(g,v0);
        ArfB0=Arf(JFF,WB0tFF,Arfdata);
	[v0,testArf]=Initv0(g,A,M,S,Frob,JFF,WB0tFF,SBasisJ2FF,WeilMatrix,Arfdata,ArfB0);print("v0 = ",lift(v0));
	B0data=[g,WeilMatrix,v0,ArfB0];
	flagArf=Checkv0(testArf,B0data);
	if(default(debug)>1,print("Check Arf: ",flagArf));
	if(flagArf==0,print("\n","ERROR: Checkv0","\n"));
	print("end");
	if(default(debug)>1,print("beginning listeOrbiteOdd"));
	[listeOrbiteOdd,listOTCNt]=InitlisteOrbiteOdd(B0data,Frob,flagTC);
	if(default(debug)>1,print("end"));
	return([listeOrbiteOdd,listOTCNt,Arfdata,B0data,Frob,flagArf]);
}

DataOTCGAdic(JInit,JFF,SBasisJ2FF,Frob,WB0tFF,WK0tGAdic,CardJq,CardJp)=
{
	my(SBasisJ2GAdic,PicTorsFrobGendata,Jfinal,WB0tGAdic);
	if(default(debug)>1,print("beginning LiftBasis"));
	\\SBasisJ2GAdic=LiftBasisJ2FF(J,SBasisJ2FF);
	[SBasisJ2GAdic,PicTorsFrobGendata]=LiftBasisJ2F2FrobGen(J,JFF,SBasisJ2FF,Frob);
	if(default(debug)>1,print("end"));
	\\print("on lift BasisJ2FF un à un: peut mieux faire");
	\\print("WE1, WE2");
	\\Jfinal=PrepChi(J,CCurve,P0,P1,P2,Frob);
	\\print("end");
	Jfinal=JInit;
	print("\n");
	print("Lift Theta Characteristic"); 
	[Jfinal,WB0tGAdic]=LiftTC(Jfinal,WB0tFF,WK0tGAdic,CardJq,CardJp);
	print("end");
	\\print("Lift de WB0tFF: Ne serez t'il pas possible de donner CardJp plutot ?");
	return([Jfinal,WB0tGAdic,SBasisJ2GAdic,PicTorsFrobGendata]);
}

InitRpzOrbiteOTC(Jf,listeOrbiteOdd,SBasisJ2GAdic,WB0tGAdic,flagTC)=
{
	my(listeRpzOrbiteOTC);
	print("\n");
	print("Computation TC under Makdisi form");
	listeRpzOrbiteOTC=InitlisteRpzOrbiteOTC(Jf,listeOrbiteOdd,SBasisJ2GAdic,WB0tGAdic,flagTC);
	print("end");
	return(listeRpzOrbiteOTC);
}
\\return([Jfinal,B0data,listeOrbiteOdd,listOTCNt,listeRpzOrbiteOTC,Arfdata,SBasisJ2GAdic,WB0tGAdic,Frob,PicTorsFrobGendata]);

/*
InitChiS(J,BasisJ2GAdic,listeOrbiteJ2,listeRpzOrbiteJ2,listeRpzOrbiteOTC)=
{
	return([InitChiJ2(J,listeRpzOrbiteJ2,listeOrbiteJ2),InitChiOTC(J,listeOrbiteJ2)]);
}
*/

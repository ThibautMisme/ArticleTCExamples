\\ \r AaTechnique.gp
\\ \r AaMakdisiMachinery.gp


InitWB0tFF(JFF,CardJp,CardJq,aMak,WK0tFF,BasisJ2,emptytors)=
\\returns the W matrix of the space representing the Point of the Jacobian B0t = B0 +(d0-(g-1))P0 - D0
{
	my(Card,char,u,v,EVEN,WB0tFF,m,N,d,flag,Weven,Wodd,uodd,vodd,dodd,projeven,projodd,projroot,WB0tFFodd,WB0tFFeven,listeRandom,G,logm,listetG,CoefWK0tFFeven,CoefWB0tFFeven);
	if(default(debug)>1,print("beginning: WB0tFF"));
	[u,v,EVEN]=gcdext(CardJp,2);
	\\print("OPTIMISATION: coder point aléatoire de J(Fp)");
	\\print("Optimisation: Gaussian elimination pour avoir J(Fq)[2^(w+1)] ou w est l'ordre de WK0tFF");
	\\print("Optimisation: transformer base J(Fp)[2^inf] en base de J(Fq)[2^inf ?");

	\\print("InitWB0tFF: BEGINNING OF THE COMPUTATIONS OVER J(FP): mettre timer","\n","THIS PART TAKES A WHILE ? TRY TO DO IT ON ANOTHER P' SUCH THAT WE CAN DO THE COMPUTATIONS ACTUALLY OVER J(FP') (but p' has to ba large enough to contain nZ points)  AND THEN SWITCH BACK OVER THE INITIAL J(FQ)","\n");
	if(default(debug)>1,print("long ?"));

	if(EVEN==1,
	\\Then there is a rational Theta Caracteristic on J(Fp)
		WB0tFF=MulInt(JFF,WK0tFF,v);
		\\B0t + (g-1)P0 is a theta caracteristiic
		);
			
	if(EVEN!=1,
		m=2^DLint(CardJp,2);\\the cardinal of the 2-Sylow of J(Fp) i.e. J(Fp)even
		N=CardJp/m;\\the cardinal of J(Fp)odd
		[u,v,d]=gcdext(N,m);
		flag=(gcd(CardJp,aMak)!=1);

		\\print("WARNING: PAS D'EXTRACTION DE RACINE SUR JP, ON PASSE A JQ");
		flag=0; \\THERE SHOULD ALWAYS BE A THETA CHARACTERISTIC ON J(Fp) SO WE SHOULD NEVER NEED TO LOOK FOR IT ON J(Fq)
		
		if( flag==0,
		\\Then there is a rational Theta Caracteristic on J(Fp), that we can find by computing the structure of G=J(Fp)[2^inf]
			Card=CardJp;
			char="Jp";		
			);
		if(flag,
		\\We have to compute B0 directly on J(Fq) by extracting a square root of WK0tFF in J(Fq)
			Card=CardJq;
			char="Jq";
			);

		m=2^DLint( Card ,2);\\cardinal of J(Fq)even or J(Fp)even
		if(default(debug)>1,print("m = ",m));
		N= Card /m;\\cardinal of J(Fq)odd OR J(Fp)odd
		if(default(debug)>1,print("N = ",N));
		[u,v,d]=gcdext(N,m);
		projeven=(u*N)% Card ;
		projodd=(v*m)% Card ;
		if(default(debug)>1,print("projeven = ",projeven));
		Weven=MulInt(JFF,WK0tFF,projeven);	
		Wodd=MulInt(JFF,WK0tFF,projodd);
			
		\\Wodd=SumJ(JFF,WK0tFF,AddFlipTM(JFF,Weven,getW0(JFF)));\\Wodd=W-Weven
                [uodd,vodd,dodd]=gcdext(N,2);
		projroot=(vodd)%N;
		WB0tFFodd=MulInt(JFF,Wodd,projroot);			
		
		if(CheckZero(JFF,Weven)!=1,
			print("\n We need to compute the even part of WK0 and therefore to produce points on J(Fp). If flagPointJFp=0, it won't work \n");
			print("emptytors = 0 = ",emptytors);	
			listeRandom=[projeven, char];
			G=StructureJFF(JFF,m,listeRandom,BasisJ2);
			logm=logint(m,2);
			listetG=InitlistetG(JFF,G);
			CoefWK0tFFeven=DLJacEven(JFF,G,listetG,logm,BasisJ2,Weven);print("CoefWK0tFFeven = ",CoefWK0tFFeven);
			CoefWB0tFFeven=Div2Liste(CoefWK0tFFeven);print("CoefWB0tFFeven =",CoefWB0tFFeven);
			WB0tFFeven=CoefBasetoWFF(JFF,CoefWB0tFFeven,FirstCol(G));
			WB0tFF=SumJ(J,WB0tFFodd,WB0tFFeven);
			);
		if(CheckZero(JFF,Weven),\\No even part
			WB0tFF=WB0tFFodd;
			);

		);

	if(CheckEqual(JFF,SumJ(JFF,WB0tFF,WB0tFF),WK0tFF)!=1,
		print("WARNING: WB0tFF wrong");
		);
	print("Test WB0tFF: ",CheckEqual(JFF,SumJ(JFF,WB0tFF,WB0tFF),WK0tFF));
	print("WB0tFF: done \n");
	return(WB0tFF);
}
	

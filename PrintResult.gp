PrintOTC(ChiOTC,flagTC)=
{
	my(a);
	print1("$\\OTC : ");
	if(flagTC,a=PrintOrbitDecomposition(factor(ChiOTC[1]),1),
		a=PrintOrbitDecomposition(factor(ChiOTC),1));
	print("$","\\","\\");
	return(a);
}

PrintETC(ChiETC,flagTC)=
{
	my(b);
	if(flagTC==0,b=0,
        	print1("$ \\ETC : ");
        	b=PrintOrbitDecomposition(factor(ChiOTC[2]),1);
        	print("$","\\","\\");
        	);
	return(b);
}

PrintJ2(ChiJ2,flagChiJ2)=
{
	my();
	if(flagChiJ2,
	        print("$ J[2] : ");
        	PrintOrbitDecomposition(factor(ChiJ2),1);
        	print("$","\\","\\");
        	);
	return();
}

PrintTHrational(a,b,listTF,listeTH,flagTC,flagTH)=
{
	my();
	print("\n");
	if((a!=0)||(b!=0), \\rational TC
        	if(a!=0,
			for(k=1,a,
				if(flagTC,print("TF rational OTC: $",polhomogenise(lift(listTF[1][k]),'z,1),"$","\\","\\"),
					print("TF rational OTC: $",polhomogenise(lift(listTF[k]),'z,1),"$","\\","\\"));
                		if(flagTH,print("Theta Hyperplane: $",lift(listeTH[k]),"$","\\","\\"));
				);
                	);
        	if(b!=0,
			for(k=1,b,
				if(flagTC,print("TF rational ETC: $",lift(listTF[2][k]),"$","\\","\\"),
					print("TF rational ETC: $",lift(listTF[2][k]),"$","\\","\\"));
                		);
			);
        	);
	return();
}

PrintFieldTC(listeFieldArchiOTC,flagTC)=
{
	my();
	print("\n");
	print("Fields of definition of the OTC","\\","\\");
	if(flagTC,for(k=1,#listeFieldArchiOTC[1],print("$",polredabs(listeFieldArchiOTC[1][k].pol),"$","\\","\\")),
		for(k=1,#listeFieldArchiOTC,print("$",polredabs(listeFieldArchiOTC[k].pol),"$","\\","\\")));
	if(flagTC,
		print("Fields of definition of the ETC","\\","\\");
		for(k=1,#listeFieldArchiOTC[2],print("$",polredabs(listeFieldArchiOTC[2][k].pol),"$","\\","\\"));
		);
	return();
}

print("Curve : $",PolCurve," = 0,$","\\","\\");
print("time computation (CPU time) : ",strtime(timecomputation),"\\","\\");
print("$P_0 = ",P0,"$,	   $P_1 = ",P1,"$,   $g = ",g,"$,   $e_{max}= ",geteJ(J),"$,   $d_0 = ",d0,"$","\\","\\");

print("$p = ",p,"$,   $q = ",p,"^",aMak,"$,   $|J(\F_p)|= ",CardJp,"$,   $|J(\F_q)|= ",CardJq,"$","\\","\\");

print("\n","Orbit Decomposition");
a=PrintOTC(ChiOTC,flagTC);
b=PrintETC(ChiETC,flagTC);

PrintTHrational(a,b,listTF,listeTH,flagTC,flagTH);

PrintFieldTC(listeFieldArchiOTC,flagTC);


Auts=[];
flagP0=1; \\if flagP0==0, we don't make a change of variable if P0 is not a pole of x
Badprime=[];
boundp=100;
P1Init=0;
P2Init=0;
flagTC=1;\\flagTC=1 if we want to compute all the TC. flagTC=0 we compute just the OTC
flagmethod="TFwithCL";\\flagmethod = "Nicolas", "Classic", "TF", "TFwithCL"
flagmethodbis="Nicolas"; \\Only if the first method failed
flagmethodTC="TFwithCL";\\flagmethod = "Nicolas", "Classic", "TF", "TFwithCL"
flagmethodTCbis="Nicolas";
flagtestbasisSymp=0; \\if flag=1 we perform an extra (time-consumming) test to check that it's okay
flagTFEonephase=1; \\flag=0 IF we don't want to try TFE in one phase. It's probably always better to try.
flagEffective=1; \\if flag=1, we try to make the TF Effective. Almost always worth the attempt.
flagmethodChiJ2="Nicolas"; \\flag= "NicolasFromScratch", "Nicolas", "Classic", 0 if we don't want ChiJ2
flagETCeffective=1;\\if flag=1 (and flagTC=1), we try to make the ETC as "effective" as we can

flagPointJFp=0;\\if flag=1, we run PicInitThibaut when there is rational torsion. If flag=0 and there is the need of computing the 2-part of WK0, the code will fail.
\\better to set flag=0 by default and change when we need it.

listPoint=[];

flagTH=1;
flagpopti=0;
flagemptytors=0;
flagETCnonVanish=0;\\if flag=1, there is a non-vanishing ETC. We perform extra computation to make the TFE better.
flagAlphaGlobal=1;\\if flagAlphaGlobal=1, Alpha is supposed defined. So we don't have to compute W2plusP0.
flagMC=0;
thresholdOrbit=150;\\above this size, we don't compute the orbit of TC
flagmethodChiJ2=0;
ChiJ2="TODO";
flagChiJ2=0;

This repository contains the code used to produce the examples of the article "Computing Theta Characteristics and their Galois Action on General Curves of Moderate Genus".<br />
Visit arXiv:2610.29449 <br />
If you enconter any issues, please contact me mismet at tcd.ie <br />



To run an example of the article:
1. Download and set up the PARI/GP language on your computed. Visit https://pari.math.u-bordeaux.fr/gp.html.<br/>
This is a necessary step as the code is implemented in PARI/GP.
2. Gather all the files of the folder 'Code' in a single folder on your computer at a certain address that I refer to as 'addressCode'.<br/>
Ex: addressCode:= ~/pari/ArticleTCExamples/Code/
2. Identify the label of the curve you are interested in. The labels are given in the Examples section of the article. They are also the name of the files in the folder 'Inputs'.<br/>
Ex: Klein
4. Open the gp interface on your computer.
5. Read the files of the folder 'Code' using the pari/gp command 'read'. You can read the file 'Code/Aaah.gp' only. It will read the other files automatically.<br/>
Ex: read("addressCode/Aaah.gp");
6. Read the file in the folder 'Inputs' that correspond to the chosen label using the command 'read'.<br/>
Ex: read("~/pari/ArticleTCExamples/Inputs/Klein.gp");
7. Read the file 'main.gp' using the command 'read'. It will run the algorithm.<br/>
Ex: read("~/pari/ArticleTCExamples/main.gp");
8. Once the algorithm has run, you can print the results as displayed in the article by reading the file 'PrintResult.gp' using the command 'read'.<br/>
Ex: read("~/pari/ArticleTCExamples/PrintResult.gp");
9. You can also directly access the data:<br/>
        - The etale polynomial(s) chi(s) are stored in the variable(s) chiOTC (and chiETC)<br/>
        - The Makdisi representations of the TC over Fq is stored as a list in the variable 'listeRpzOrbiteOTC'.<br/>
            This list is split in two if the ETC have been computed. It is furthermore split by Qq-orbits. Only one representative by orbit is stored.<br/>
        - The Theta Fractions are stored as a list in the variable 'listTF'.<br/>
            This list is split in two if the ETC have been computed. It is furthermore split by Q-orbits. Only one representative by orbit is stored.<br/>
        - The Theta Hyperplanes are stored as a in the variable 'listTH'.<br/>
            This list is plit by Q-orbits. Only one representative by orbit is stored.<br/>
        - The fields of definition of the TC are stored as a list in the variable 'listFieldTC'. They are as returned by the algorithm. In particular, the commend 'polredabs' has not been applied to them.<br/>
            This list is split in two if the ETC have been computed.


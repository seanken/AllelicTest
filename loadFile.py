######################
##Loads Plink file
######################
import numpy as np;
from pysnptools.snpreader import Bed;
from pysnptools.snpreader import Pheno;
from os.path import isfile;

##
##If D is a n by m matrix returns m by 3 array r,
##where r[i][j] are the number of entries in column i of D equal to j
##
def getMarginals(D):
	n=len(D);
	tot=np.sum(D,axis=0);
	tot2=np.sum(np.minimum(D,1),axis=0);
	r2=tot-tot2;
	r1=tot2-r2;
	r0=n-tot2;
	return [[r0[i],r1[i],r2[i]] for i in range(0,len(r0))];

##
##gets y, the phenotype, X, the normalized matrix, sFil,
##the snpreader, and Q the covariate matrix
##
def getData(filename):
	mph=3;
	sFil=Bed(filename);
	yFil=Pheno(filename+".fam");
	snpList=sFil.sid;
	y=yFil.read().val[:,mph];
	y=[i-1 for i in y]
	Icases=[i for i in range(0,len(y)) if y[i]>0];
	Icont=[i for i in range(0,len(y)) if y[i]<1];
	sFilcases=sFil[Icases,:]
	sFilcont=sFil[Icont,:]


	Dcont=sFilcont.read().val;
	Dcases=sFilcases.read().val;

	
	r=getMarginals(Dcont);
	s=getMarginals(Dcases);

	return [r,s,snpList];
	

if __name__=="__main__":
	D=[[0,1,0],[0,2,0],[0,0,0],[1,1,1]]
	print getMarginals(D)

	filename="../PrivGWAS/GWAS/cleaned";
	[r,s,slist]=getData(filename);
	print slist[:10];
	print len(r);
	print len(s);
	print sum(r[0]);
	print sum(s[0]);

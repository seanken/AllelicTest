##############################################
##This implements the neighbor distance method for picking high scoring SNPs in a 
##privacy preserving manner. Unlike previous approaches this takes constant time per SNP.
##Note this is still a fairly early version.
################################################
##Notation:
##R number of cases
##S number of controls
##N=R+S
##c is threshold
##r0,r1,r2 number of cases with 0, 1 or 2 copies of minor allele
##s0,s1,s2 same for controls
##x=2r0+r1
##y=2s0+s1

##
##Dependencies
##
import math;
import numpy as np;
import scipy as sp;
import random;
from numpy.random import multinomial as multi;
from numpy.random import laplace as Lap;

##
##Calc set Q (boundary points in optimization)
##
def getQ(r0,r1,r2,s0,s1,s2,R,S,N,c):
	Q=[];
	##Want to find all (p0,p1) that should be in Q
	for p0 in [2*(r0+r2)+r1,2*r0+r1,r1,0,2*R]:##first consider fixed p0
		p=float(p0);
		c0=2*N/float(c*R*S);
		c0=1/c0;
		##use quadratic equation to find possible values for p1
		A=R**2/c0+1;
		B=2.0*p-2*N-2*p*R*S/c0;
		C=p**2*S**2/c0+p**2-2*N*p;
		det=B**2-4*A*C;
		if det<0:
			continue;
		det=math.sqrt(det);
		p1p=(-B+det)/float(2*A)	
	
		p1n=(-B-det)/float(2*A)	
		if p1p<=2*S and p1p>=0 and p+p1p!=0 and p1p+p!=2*N:
			Q.append([p,p1p]);

		if p1n<=2*S and p1n>=0 and p+p1n!=0 and p1n+p!=2*N:
			Q.append([p,p1n]);

	for p1 in [2*(s0+s2)+s1,2*s0+s1,s1,0,2*S]:##then similar for fixed p1
		p=float(p1);	
		c0=2*N/float(c*R*S);
		c0=1/c0;
		A=S**2/c0+1;
		B=2.0*p-2*N-2*p*R*S/c0;
		C=p**2*R**2/c0+p**2-2*N*p;
		det=B**2-4*A*C;
		if det<0:
			continue;
		det=math.sqrt(det);
		p1p=(-B+det)/float(2*A)	
	
		p1n=(-B-det)/float(2*A)	

		if p1p<=2*R and p1p>=0 and p1p+p!=0 and p1p+p!=2*N:
			Q.append([p1p,p]);

		if p1n<=2*R and p1n>=0 and p1n+p!=0 and p1n+p!=2*N:
			Q.append([p1n,p]);

	return Q;


##
##Calc set P, points with derivatives in 0,1,2,.5
##
def getP(R,S,N,c):
	P=[];
	for dydx in [.5,2.0,1.0,0.0]:##iterates through all possible derivatives
		c0=c*R*S/float(2*N);
		a=1-R*S/c0+(1+(R*R)/c0)*dydx;
		b=1+(S*S)/c0+(1-R*S/c0)*dydx;
		d=-N*(dydx+1);
		u=-b/a;
		v=-d/a;
		##use quad equation to find possible points
		A=S*S-2*u*R*S+u*u*R*R+c0*(1+u)**2;
		B=-2*v*R*S+2*u*v*R**2+c0*(1+u)*(v-2*N)+c0*v*(1+u);
		C=v**2*(R*R)+c0*v*(v-2*N);
		B=B/A;
		C=C/A;
		A=1.0;
		
		det=B*B-4*A*C;
		det=math.sqrt(det);
		x1=(-B+det)/float(2*A);
		y1=u*x1+v
		if 0<=x1 and x1<=2*R and 0<=y1 and y1<=2*S:
			P.append([x1,y1]);##adds first one to P
		

		x2=(-B-det)/float(2*A)

		y2=u*x2+v
		if 0<=x2 and x2<=2*R and 0<=y2 and y2<=2*S:
			P.append([x2,y2]);##adds second one to P
	return P;
		
##
##Helper function to calc g1 or g2
##
def calcGi(x,r0,r1,r2,R):
	if 2*(r0+r2)+r1>=x and x>= 2*r0+r1:
		return (x-2*r0-r1)/2.0;
	if r1<=x and x<=2*r0+r1:
		return (2*r0+r1-x)/2.0;
	if x>=2*(r0+r2)+r1:
		return r2+x-2*(r0+r2)-r1;
	return r0+r1-x;

##
##calculates the relaxed function g, that approximates the function we want to optimize
##
def calcG(x,y,r0,r1,r2,s0,s1,s2,R,S):
	return calcGi(x,r0,r1,r2,R)+calcGi(y,s0,s1,s2,S);	
		


##
##Checks if database at distance delta from current database with score <=c.
##
def CheckDelta(delta,r0,r1,r2,s0,s1,s2,R,S,N,c):
	x=2*r0+r1;
	y=2*s0+s1;
	##Sanity Check
	if 2*R<x or x<0 or y<0 or 2*S<y:
		return False;
	
	##Sets describing affine maps that need to iterate over
	U=[[2*r0+r1,2*(r0+r2)+r1],[r1,2*r0+r1],[2*(r0+r2)+r1,2*R],[0,r1]];
	U=[U[0],U[0],U[1],U[1],U[2],U[2],U[3],U[3]]
	V=[.5,.5,-.5,-.5,1.0,1.0,-1.0,-1.0];
	D=[(-2*r0-r1)/2.0,(-2*r0-r1+1)/2.0,(2*r0+r1)/2.0,(2*r0+r1+1)/2.0,-r2-2*r0-r1,-r2-2*r0-r1,r0+r1,r0+r1];
	C=[(-2*r0-r1)%2,(-2*r0-r1+1)%2,(-2*r0-r1)%2,(-2*r0-r1+1)%2,1,0,1,0];
	Up=[[2*s0+s1,2*(s0+s2)+s1],[s1,2*s0+s1],[2*(s0+s2)+s1,2*S],[0,s1]];
	Up=[Up[0],Up[0],Up[1],Up[1],Up[2],Up[2],Up[3],Up[3]]
	Vp=[.5,.5,-.5,-.5,1.0,1.0,-1.0,-1.0];
	Dp=[(-2*s0-s1)/2.0,(-2*s0-s1+1)/2.0,(2*s0+s1)/2.0,(2*s0+s1+1)/2.0,-s2-2*s0-s1,-s2-2*s0-s1,s0+s1,s0+s1];
	Cp=[(-2*s0-s1)%2,(-2*s0-s1+1)%2,(-2*s0-s1)%2,(-2*s0-s1+1)%2,1,0,1,0];


	## iterate over all pairs
	for i in range(0,len(U)):
		for j in range(0,len(Up)):
			##parameterization of lines
			alpha=[1.0,-V[i]/Vp[j]];
			beta=[0,(delta-D[i]-Dp[j])/Vp[j]];
			
			t1=[(U[i][0]-beta[0])/alpha[0],(U[i][1]-beta[0])/alpha[0]];
			t2=[(Up[j][0]-beta[1])/alpha[1],(Up[j][1]-beta[1])/alpha[1]];
			if t2[0]>t2[1]:
				t2=[t2[1],t2[0]]
			if t1[0]>t1[1]:
				t1=[t1[1],t1[0]]
			a=max(t1[0],t2[0]);
			b=min(t1[1],t2[1]);
			if b<a:
				continue;
			##use quadratic equation to find interval so parameterization gives solution in that interval
			c0=c*R*S/float(2*N);
			A=alpha[0]**2*S**2+alpha[1]**2*R**2-2*alpha[0]*alpha[1]*R*S+c0*(alpha[0]+alpha[1])**2;
			B=-2*alpha[0]*beta[1]*R*S+2*alpha[1]*beta[1]*R**2+c0*(alpha[0]+alpha[1])*(2*beta[1]-2*N);
			Ctemp=beta[1]**2*R**2+c0*beta[1]*(beta[1]-2*N);
			A=float(A);
			B=float(B);
			Ctemp=float(Ctemp);
			if A<0:
				A=-A;
				B=-B;
				Ctemp=-Ctemp;
			det=B**2-4*A*Ctemp;
			if det<0:
				continue;
			det=math.sqrt(det);
			p0=(-B-det)/float(2*A)
			x0=alpha[0]*p0+beta[0];
			y0=alpha[1]*p0+beta[1];
			p1=(-B+det)/float(2*A)
			a=max(a,p0);
			b=min(b,p1);
			if b<a:	
				continue;
			if math.ceil(a)>b:
				continue;

			##Check if any integer solution in interval meeting modularity constraints-- if so return it, else continue
			ca=int(math.ceil(a));
			if not (alpha[1]*ca+beta[1]).is_integer():
				ca=ca+1;
			if ca>b:
				continue;
			if int(alpha[0]*ca+beta[0])%2==C[i] and int(alpha[1]*ca+beta[1])%2==Cp[j]:
				return [alpha[0]*ca,alpha[1]*ca+beta[1]];
			ca=ca+1;
			if not (alpha[1]*ca+beta[1]).is_integer():
				ca=ca+1;
			if ca>b:
				continue;
			if int(alpha[0]*ca+beta[0])%2==C[i] and int(alpha[1]*ca+beta[1])%2==Cp[j]:
				return [alpha[0]*ca,alpha[1]*ca+beta[1]];

	return [];##if no such points exist








##
##Calculates the neighbor distance using our algorithm
##
def neighDistSig(r0,r1,r2,s0,s1,s2,c,R,S,N,P=[],rnd=False):
	x=2*r0+r1;
	y=2*s0+s1;
	##Some sanity checks
	if c<=0:
		return -1;
	if r0<0 or r1<0 or r2<0 or s0<0 or s2<0 or s1<0:
		return -1;
	if len(P)==0:##if did not already generate P
		P=getP(R,S,N,c);

	Q=getQ(r0,r1,r2,s0,s1,s2,R,S,N,c);##Calculates set Q

	delta1=min([calcG(p[0],p[1],r0,r1,r2,s0,s1,s2,R,S) for p in P])##min distance in P

	delta2=min([calcG(p[0],p[1],r0,r1,r2,s0,s1,s2,R,S) for p in Q])##min distance in Q
	delta=min(delta1,delta2);
	delta=math.ceil(delta);##optimum of relaxed problem
	if rnd:##called if not a significant SNP
		return delta;

	##now round to integer solution
	for d in [delta-2,delta-1,delta,delta+1,delta+2,delta+3]:
		keyDel=CheckDelta(d,r0,r1,r2,s0,s1,s2,R,S,N,c)##Checks if exists database within distance d with score less than c
		if len(keyDel)>0:##if such a database exists
			return d;

	return -1;##Means failed








##
##Calculates neighbor distance
##
def neighDist(r0,r1,r2,s0,s1,s2,c,R,S,N,P,simp=False):
	x=2*r0+r1;
	y=2*s0+s1;
	if c<=2*N*(x*S-y*R)**2/float(R*S*(x+y)*(2*N-x-y)):
		dst=neighDistSig(r0,r1,r2,s0,s1,s2,c,R,S,N,P);##if greater than boundary	
		return dst;
	dst=neighDistSig(r0,r1,r2,s0,s1,s2,c,R,S,N,rnd=True);##if less than boundary
	if simp:
		return dst;
	return 1-dst;


##
##Calculates neighbor distance for all SNPs
##With the given boundary c
##r is an ordered list of 3 tuples, one for each SNP, corresponding to r0,r1 and r2
##Similarly for s, except with s0,s1,s2
def NeighDistTot(r,s,c,simp=False):
	m=len(r);
	R=r[0][0]+r[0][1]+r[0][2];##number of cases
	S=s[0][0]+s[0][1]+s[0][2];##number of controls
	N=R+S;##number of participants
	P=getP(R,S,N,c);##Set P used in calculation
	return [neighDist(r[i][0],r[i][1],r[i][2],s[i][0],s[i][1],s[i][2],c,R,S,N,P,simp=simp) for i in range(0,m)];




##
##Calculates the allelic test statistic
##For each SNP
##
def AllelicScoreTot(r,s):
	m=len(r);
	R=r[0][0]+r[0][1]+r[0][2];##number of cases
	S=s[0][0]+s[0][1]+s[0][2];##number of controls
	N=R+S;##number of participants
	return [2*N*((2*r[i][0]+r[i][1])*S-(2*s[i][0]+s[i][1])*R)**2/float(R*S*(2*r[i][0]+r[i][1]+2*s[i][0]+s[i][1])*(2*N-2*r[i][0]-r[i][1]-2*s[i][0]-s[i][1])) for i in range(0,m)];


##
##Picks mret SNPS using the Laplacian bases approached
##using given scores, epsilon and sensitivity
##
def LapPick(mret,epsilon,scores,sens):
	m=len(scores)
	scoresLap=[s+Lap(0.0,sens*mret*2/epsilon) for s in scores];##perturbed scores
	return sorted([i for i in range(0,m)],key=lambda i:-scoresLap[i])[:mret];



##
##Picks mret SNPS using the exponential mechanism based approached
##using given scores, epsilon and sensitivity
##
def ExpPick(mret,epsilon,scores,sens):
	m=len(scores);
	sc=scores[0];
	neighs=[epsilon*s/(mret*2.0*sens) for s in scores];
	mxn=max(neighs);
	minNeigh=min(neighs);
	neighs=[nei-mxn for nei in neighs];
	ret=[];##To return
	##picks mret SNPs, one by one, with no repetition

	neighs=[math.exp(nei) for nei in neighs];
	Z=sum(neighs);
	neighs=[nei/Z for nei in neighs];##probabilities for each SNP
	for i in range(0,mret):
		mlt=multi(1,neighs)##picks next SNP to return
		k=min([j for j in range(0,m) if mlt[j]>0]);
		scores[k]=min(scores);
		ret.append(k);
		##have to recalc probs to avoid underflow
		con=epsilon/(mret*sens*2.0)
		neighs=[con*float(s) for s in scores];
		mxn=max(neighs);
		neighs=[nei-mxn for nei in neighs];
		neighs=[math.exp(nei) for nei in neighs];
		for k in ret:
			neighs[k]=0.0;
		Z=sum(neighs)
		neighs=[nei/Z for nei in neighs]	
	return ret;




##
##Perturb the input, helper method for calcBoundary
##
def pertData(r0,r1,r2,s0,s1,s2,epsilon=1.0):
	R=r0+r1+r2;
	S=s0+s1+s2;
	N=R+S;
	x=2*r0+r1;
	y=2*s0+s1;
	xdp=x+Lap(0.0,2/epsilon);##perturbed x
	ydp=y+Lap(0.0,2/epsilon);##perturbed y
	return max(2*N*((xdp)*S-(ydp)*R)**2/float(R*S*(ydp+xdp)*(2*N-xdp-ydp)),2*N/float(2*N-1));

##
##Calculates boundary for neighbor distance
##
def calcBoundary(r,s,mret,epsilon):
	scores=AllelicScoreTot(r,s);
	##finds mret and mret+1st highest scoring SNPs
	atBoundary=sorted([i for i in range(0,len(r))], key=lambda i:-scores[i])[mret-1:mret+1];
	i1=atBoundary[0];
	i2=atBoundary[1];
	threshDP=pertData(r[i1][0],r[i1][1],r[i1][2],s[i1][0],s[i1][1],s[i1][2],epsilon=epsilon/2.0)+pertData(r[i2][0],r[i2][1],r[i2][2],s[i2][0],s[i2][1],s[i2][2],epsilon=epsilon/2.0);
	return threshDP/2.0


##
##Applies exponential mechanism using neighbor score
##
def neighExp(mret,epsilon,r,s):
	sens=1.0;
	c=calcBoundary(r,s,mret,.1*epsilon);
	scores=NeighDistTot(r,s,c)
	return ExpPick(mret,.9*epsilon,scores,sens)


##
##Applies exponential mechanism using neighbor score
##with arbitrary c
##
def neighExpArb(mret,epsilon,r,s,c):
	sens=1.0;
	e=epsilon;
	scores=NeighDistTot(r,s,c)
	return ExpPick(mret,epsilon,scores,sens)

##
##Applies Laplace Method using neighbor score
##
def neighLap(mret,epsilon,r,s):
	sens=1.0;
	c=calcBoundary(r,s,mret,.1*epsilon);

	scores=NeighDistTot(r,s,c)
	return LapPick(mret,.9*epsilon,scores,sens)

##
##Applies exponential mechanism using naive score
##
def scoreExp(mret,epsilon,r,s):
	R=r[0][0]+r[0][1]+r[0][2];
	S=s[0][0]+s[0][1]+s[0][2];
	N=R+S;
	scores=AllelicScoreTot(r,s);
	##calc sensitivity (See Uhlers paper)
	sens1=8*N*N*S/float(R*(2*S+3)*(2*S+1));
	sens2=4*N**2*((2*R**2-1)*(2*S-1)-1)/float(R*S*(2*R+1)*(2*R-1)*(2*S+1));
	sens3=8*N**2*R/float(S*(2*R+3)*(2*R+1));
	sens4=4*N**2*((2*S**2-1)*(2*R-1)-1)/float(R*S*(2*S+1)*(2*S-1)*(2*R+1));
	sens=max(sens1,sens2,sens3,sens4);
	
	return ExpPick(mret,epsilon,scores,sens);

##
##Applies Laplacian based method using naive score
##
def scoreLap(mret,epsilon,r,s):	
	R=r[0][0]+r[0][1]+r[0][2];
	S=s[0][0]+s[0][1]+s[0][2];
	N=R+S;
	scores=AllelicScoreTot(r,s);

	##Calc sensitivity
	sens1=8*N*N*S/float(R*(2*S+3)*(2*S+1));
	sens2=4*N**2*((2*R**2-1)*(2*S-1)-1)/float(R*S*(2*R+1)*(2*R-1)*(2*S+1));
	sens3=8*N**2*R/float(S*(2*R+3)*(2*R+1));
	sens4=4*N**2*((2*S**2-1)*(2*R-1)-1)/float(R*S*(2*S+1)*(2*S-1)*(2*R+1));
	sens=max(sens1,sens2,sens3,sens4);

	return LapPick(mret,epsilon,scores,sens);


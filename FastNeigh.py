###
##User interface
#########
import PrivPick;
import sys;
import loadFile;

def UI(filename,mret,eps):
	print "Lets pick some SNPs!"
	print "The number of SNPs to return is: "+str(mret);
	print "The file to look at is: "+filename;
	print "The value of epsilon is: "+str(eps);
	print "\n";
	print "Getting data"
	[r,s,SNPList]=loadFile.getData(filename);
	print "Got data!\n\n"
	print "Picking SNPs"
	I=PrivPick.neighExp(mret,eps,r,s);
	print "Got SNPs!"
	print "The SNPs are: "
	for i in I:
		print SNPList[i];
	print "Good Bye!"


if __name__=="__main__":
	args=sys.argv;
	if len(args)!=4:
		print "Not enough arguments!"
	else:
		mret=int(args[2]);
		eps=float(args[3]);
		filename=args[1];
		UI(filename,mret,eps);

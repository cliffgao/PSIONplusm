#!/usr/bin/env python
"""
to compute the confuse matrix and auc
inputfile: >probability; number; modify 2012-10-18
change in Binary type 2015-05-26
"""
import os,math,sys
#import os,math,sys,numpy,scipy

###	
def  get_prob_based_label(probs,k):
	N=len(probs)
	results=[]
	for i in range(N):
		results.append(probs[i][k])
	return results
	
###
def main(fr_n,Binary,cutoff=0.5):  
  # store the results
	labels=[] #actual labels
	probs=[] #probablities
	#read in svm-output files
	fr=open(fr_n,'r')
	lines=[line.strip() for line in fr.readlines()]
	fr.close() # 1st id; 2nd probability 3rd actuals
	#probs=[]
	MCCS=[]
	if Binary=="T":
		for i in range(len(lines)):
			if len(lines[i])>0:
				aline=lines[i].split()
				alabel=aline[0]
				aprob=aline[-1] # order -1 1; -1 is potive label
				probs.append(float(aprob))
				labels.append(alabel)
		P_label="1"		
		##MCCS=grid_mcc(labels,probs,P_label)		
		mcc,acc,fm,tp,fn,fp,tn,Sn,Sp=get_mcc(labels,probs,cutoff,P_label)
		print("%s cutoff|fmeasure|mcc|acc|sn|sp" %(fr_n))
		print("%.3f: %.3f,%.3f,%.3f,%.3f,%.3f" %(cutoff,fm,mcc,acc,Sn,Sp))
	else:#multi-class
		for i in range(len(lines)):
			if len(lines[i])>0:
				aline=lines[i].split()
				alabel=aline[0]
				# store 1 2 3 4 probs
				aprobs=[]
				for j in range(1,len(aline)):
					aprobs.append(float(aline[j]))
					
				probs.append(aprobs)
				labels.append(alabel)
	
		P_labels=["1","2","3","4" ]
		#MCCS=[]
		tmpmcc=[]
		tmpacc=[]
		tmpfm=[]
		tmptp=[]
		for k in range(len(P_labels)):
			P_label=P_labels[k]
			pprob=get_prob_based_label(probs,k)
			#sp1,sn=roc_xy_axis(labels,results,cutoff,P_label,N_label)
			xmcc,xacc,xfm,xtp,xfn,xfp,xtn,xSn,xSp=get_mcc(labels,pprob,cutoff,P_label)
			tmpmcc.append(xmcc)
			tmpacc.append(xacc)
			tmpfm.append(xfm)
			tmptp.append(xtp)
			print("Label: %s %.3f" %(P_label, cutoff))
			print("%.3f,%.3f,%.3f,%.3f,%.3f" %(xfm,xmcc,xacc,xSn,xSp))
			#mccsum=sum(tmpmcc)
			#if mccsum=0
		mcc=sum(tmpmcc)/len(tmpmcc)
		acc=sum(tmpacc)/len(tmpacc)
		fm=sum(tmpfm)/len(tmpfm)
		oacc=sum(tmptp)*1.0/len(labels)
		print("AVERAGE: %s cutoff|fmeasure|mcc|acc|overall" %(fr_n))
		print("AVERAGE:%.3f: %.3f,%.3f,%.3f,%.3f" %(cutoff,fm,mcc,acc,oacc))
		
			#print cutoff  #2010-09-01
	

	
	
"""	
	fw_n=fr_n+".auc.confuse.matrix.txt"
	fw=open(fw_n,'w')
	fw.write("#labels predictions\n")
	if len(results)!=len(labels):
		print ("%s results and labels length is different")
	for i in range(len(results)):
		fw.write("%s %9f\n" %(labels[i],results[i]))

	xyaxis=[]
	#print xaxis
	xyaxis.sort()
	#print xyaxis #2010-09-01
	fw.write("#xyaxis ==ROC x-axis:1-Sp y-axis:Sn==\n")  
	for k in range(len(xyaxis)):
		xa,yb=xyaxis[k]
		fw.write("%f   %f\n" %(xa,yb))
	myauc=compute_auc(xyaxis)
	print "AUC= ", myauc
	fw.write("AUC=%.3f\n"  %myauc)
	#compute the confuse matrix
	cutoff=0.5
	(tp,fn,fp,tn,p,n,ACC,Sn,Sp,MCC,Precision,F_measure)=roc_xy_axis(labels,results,cutoff,P_label,N_label)
	#write the confuse matrix
	fw.write("#matrix ------Confusion  Matrix-----\n")
	fw.write("TP FP|%5d %5d\n" %(tp,fp)) 
	fw.write("FN TN|%5d %5d\n" %(fn,tn)) 
	fw.write("Accuracy            %.3f \n" %ACC)
	fw.write("Sn/TPR/recall       %.3f \n" %Sn)
	fw.write("Sp/TNR              %.3f \n" %Sp)
	fw.write("MCC                 %.3f \n" %MCC)
	fw.write("Precision/PPR       %.3f \n" %Precision)
	fw.write("F_measure           %.3f \n" %F_measure)
	fw.close()

"""
#fr_n="cv5.cdhit30.Types.vote.1.prob.txt"
#Binary="F"
#main(fr_n,Binary)
if __name__=="__main__":
  #comput_svm(float(sys.argv[1]),float(sys.argv[2]))
	if len(sys.argv)<2:
		print("Usage: my_auc_confuse.py  fr_n Binary(T/F) cutoff")
	else:
		fr_n=sys.argv[1]
		Binary=sys.argv[2]
		cutoff=float(sys.argv[3])
		main(fr_n,Binary,cutoff)
  

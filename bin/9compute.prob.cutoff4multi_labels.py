#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 19:19:23 2015
#2015-05-01  compute just for acc =tp+tn / N
@author: cliff
"""
import os,sys
EPS=0.0000001
# 2015-04-11 add the predictions from PSION
#import numpy as np 
def get_orders(oneline):
	header=oneline.split()
	orders=[]
	#if len(header) >3 : #multi-classes
	if header[0]=="labels":
		if len(header) >3:
		#multi-classes  1 2 3 4
			for i in range(4):
				for j in range(1,len(header)):
					if int(header[j])== (i+1):
						orders.append(j)
		elif len(header)==3:
		# binary-classes label  -1  1  # begain with label -1
			if header[1]=="1":
				orders.append(2)
				orders.append(1)
			elif header[2]=="1": # label -1 1
				orders.append(1)
				orders.append(2)

	else:
		print("not right head")
	return orders


def get_predict_from_psion(fr_n):
	#read in svm
	fr=open(fr_n,'r')
	lines=[line.strip() for line in fr.readlines()]
	fr.close()
	
	labels=[]
	predicts=[]
	firstLine=lines[0]
	orders=get_orders(firstLine)
	for k in range(1,len(lines)):
		if len(lines[k])>0:
			aprob=[]
			aline=lines[k].split()
			alabel=aline[0]
			for idx in orders:
				aprob.append(float(aline[idx]))
			#for j in range(1,len(aline)):
			#	aprob.append(float(aline[j]))
			labels.append(alabel)
			predicts.append(aprob)
	return (labels,predicts)
			

	

def add2lists(alist,blist):
	clist=[0.0 for k in range(len(alist))]
	if len(alist)!=len(blist):
		print("no same length")
		print(alist)
		print(blist)
	else:
		for k in range (len(alist)):
			clist[k]=(alist[k]+blist[k])*1.0/2
	return clist
def merge_psion_blast(svm_label,qidx,label,svm_predict,MAP,counterP,counter):
	#counterP: probability from blast
	#if blast return noting using psion prediction
	probs=[]
	N=len(counterP)
	probs=[0.0 for k in range(N)]
	#ncount=0
	if svm_label[qidx]!=label:
		#print("PSION label: %s blast label:%s" %(svm_label[qidx],str(label)))
		pass
	if True : #label is same , combine counter
		#print(qidx)
		asvmP=svm_predict[qidx]
		NoHit=False
		for k,v in counter.items():
			if k=='N'and v==1:
				NoHit=True
		if NoHit: # if no hit using psion prediction
			probs=asvmP
			#ncount=ncount+1
			#print("no hit", probs)
		else: #else using BLAST prediction
			#probs=add2lists(asvmP,counterP)
			probs=counterP
			#print("hit",probs)
		#for k in range(N):
		#	probs[k]=(counterP[k]+asvmP[k])/2.0
	#if ncount>0:
	#	print(ncount)
	return (probs)

	
def prob2label(probs,Binary):
	#labels=["1","2","3","4" ]
	pmax=max(probs)
	try:
		idx=probs.index(pmax)
	except:
		print("Can not find max prob")
		idx=-1
	if idx>=0:
		if Binary=="T": #begin -1 1
			if idx==0:
				label="-1"
			elif idx==1:
				label="1"
			else:
				print("no such label")
		else: #multi-class  1 2 3 4
			label=str(idx+1)
		
	else:
		print("max prob not find")
		label="5"
	return label	
def prob2label4blast(counter,pcounter,Binary):
	#labels=["1","2","3","4" ]
	#print("counter",counter)
	#print("pcounter",pcounter)
	for k,v in counter.items():
		if k=="lig" or k=="vol":
			model="VLG"
			break
		if k=="ion" or k=="non_ion":
			model="ION" #voltage-gated subtypes
			break
		if k=="k" or k=="anion" or k=="ca" or k=="na":
			model="VGS"#voltage-gated subtypes
			break
		
	if Binary=="T":
		probs=[0,0] #2 labels
	else:
		probs=[0,0,0,0]# 4label
	if counter['N']>0: #nohit
		if Binary=="T":
			if model=="ION":
				probs=[1.0,0.0] # nohit is non-ionchannel (-1)
			elif model=="VLG":
#				probs=[0.0,1.0] # nohit is voltage-gated
				probs=[1.0,0.0] # nohit is ligand-gated (-1)
				print("model VLG")
		else:
			probs=[1.0,0.0,0.0,0.0]
	else:
		for k,v in counter.items():
			if v==0:
				pcounter[k]=0.0
			else:
				pcounter[k]=pcounter[k]*1.0/counter[k]
			if k=="lig" or k=="non_ion" or k=="k":# negative 
				probs[0]=pcounter[k]
			if k=="vol" or k=="ion" or k=="anion":#positive
				probs[1]=pcounter[k]
			if k=="ca":
				probs[2]=pcounter[k]
			if k=="na":
				probs[3]=pcounter[k]				
	#print("pcounter",pcounter)	
	#normalized the probs print("pcounter",pcounter)	
	xtmp=sum(probs)#probs=np.average(probs)#normalized the probs print("pcounter",pcounter)	
	if abs(xtmp)<EPS:
		pass
	else:
		for j in range(len(probs)):
			probs[j]=probs[j]*1.0/xtmp
		
	pmax=max(probs)
	try:
		idx=probs.index(pmax)
	except:
		print("Can not find max prob")
		idx=-1
	if idx>=0:
		if Binary=="T":
			if idx==0:
				label="-1"
			elif idx==1:
				label="1"
			else:
				print("no such label")
		else: #multi-class
			label=str(idx+1)
		
	else:
		print("max prob not find")
		label="5"
	return (label,probs)		
def get_label_prediction(fr_n,fw_n,MAP,psion_fn,Ecutoff,blastOutFlag=False):
	# get the labe from fast files >id;blast1,blast2
	fr=open(fr_n,'r')
	lines=[line.strip() for line in fr.readlines()]
	fr.close()
	if len(MAP)==3:
		Binary="T"
	else:
		Binary="F"
	# counter the vote
	counter={}
	#for k,v in MAP.items():
	#	counter[k]=0
	#---read in the psion-(svm model )predictions
	svm_label,svm_predict=get_predict_from_psion(psion_fn)
	fw=open(fw_n,'w') # combined psion_results and results
	fw.write("#ION: ion (1) non-ion(-1);VLG: ligand-gated(-1) voltage-gated(1); VGS: potassium(K,1) anion(Anion 2) calcium(Ca 3) sodium(Na 4)\n")
	if blastOutFlag:
		fw2=open(fw_n+".justblast",'w') #write the blast prediction
		fw2.write("#ION: ion (1) non-ion(-1);VLG: ligand-gated(-1) voltage-gated(1); VGS: potassium(K,1) anion(Anion 2) calcium(Ca 3) sodium(Na 4)\n")
	ncount=0
	qidx=0
	for i in range(len(lines)):
		if len(lines[i])>0 and lines[i][0]=='>':
			firstline=lines[i][1:].split("|")
			proid=firstline[0]
			#label=MAP[proid]
			label="-1"
			#seq=lines[i+1]
			counter={} #count how many hits from blast 
			pcounter={} #accumulate probs
			for k,v in MAP.items():
				counter[k]=0
				pcounter[k]=0.0

			j=i+1
			seq=lines[j]
			if seq[0]=='-': #blast No hit
				predict='N'
				BlastScore=1
				ncount=ncount+1 # sum the no-hit number
			else:
				predict=seq.split("|")[0]
				Evalue=float(seq.split()[-1])
				BlastScore=Ecutoff*1.0/(Ecutoff+Evalue)
			counter[predict]=counter[predict]+1
			pcounter[predict]=pcounter[predict]+BlastScore


			#blast prediction	
			apredict,blastProb=prob2label4blast(counter,pcounter,Binary)
			if blastOutFlag:
				fw2.write("%s\n" %(apredict))
			#if predict=="N":
			#	print("blast:",apredict)
			#combine -PSION results
#			bp: blast +psion
			label="-1" 
			bpProb=merge_psion_blast(svm_label,qidx,label,svm_predict,MAP,blastProb,counter)
			#print(blastProb)
			#print(svm_label[qidx],svm_predict[qidx])
			#print(bpProb)
			xapredict=prob2label(bpProb,Binary)#psion_plus_blast
			fw.write("%s " %(xapredict))
			for xitem in bpProb:
				fw.write("%f " %xitem)
			fw.write("\n")
			#print(blastProb,bpProb,xapredict)
			#if predict=="N":
			#	print("Line: %d label: %2s, blast: %s svm: %s" %(qidx+1,label,apredict,xapredict))
			#if qidx==59:
			#	print(label,apredict,xapredict)


			
			qidx=qidx+1
	#print("ecutoff,No-hit number:", Ecutoff,ncount)		
	fw.close()
	if blastOutFlag:
		fw2.close()
	return(0)
	


if __name__=="__main__":
	if len(sys.argv)<3:
		print("Usage get_label_prediction() fr_n(blast fasta) fw_n (label/predict) model(Ion_Non,Vol_Lig,Types) psion_out")
	else:
		fr_n=sys.argv[1] #"myblast.Vol_Lig.multi.fa"
		fw_n=sys.argv[2] #"eg.myblast.Vol_Lig.label.predict.txt"
		model=sys.argv[3]
		psion_fn=sys.argv[4] #"out.vol_lig.txt"
		Ecutoff=float(sys.argv[5])
		blastOutFlag=sys.argv[6]
		if blastOutFlag=="True":
			blastOutFlag=True
		else:
			blastOutFlag=False

		if model=="1": #"Ion_Non":
			MAP={'non_ion':'-1','ion':'1','N':'2'}
			Binary="T"
		elif model=="2": #"Vol_Lig":
			MAP={'lig':'-1','vol':'1','N':'2'}
			Binary="T"
		elif model=="3": #"Types":
			MAP={'k':'1','anion':'2','ca':'3','na':'4','N':'5'}
			Binary="F"
		else:
			print("model should be 	Ion_Non,Vol_Lig,Types")
			MAP=None
		if MAP is None:
			print("Wrong model type")
		else:
			get_label_prediction(fr_n,fw_n,MAP,psion_fn,Ecutoff,blastOutFlag)


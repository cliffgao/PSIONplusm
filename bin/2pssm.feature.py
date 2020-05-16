#!/usr/bin/env python
"""
this scripe is uesd to compute the pssm features 
2010-12-20 
"""
#from my_add_two_vectors import My_add_two_vectors
#from my_pssm_normalize import My_pssm_normalize
import os
import sys
import math

"""
normailize the pssm socre
using 
1
-------
1+exp(-x)
"""
def My_pssm_normalize(myvector):
	results=[]
	for i in range(len(myvector)):
		ori_item=float(myvector[i])
		try:
			new_item=1/(1+math.exp(-1*ori_item))
		except OverflowError:
			new_item=0
			#print"#",ori_item,"#",myvector
		results.append(new_item)
	return (results)


def r_array(file_name):
	"""read in one pssm file return an array and a fasta sequence"""
	myarray=[]
	myseq=""
	try:
		myfile=open(file_name,'r')
	except IOError:
		print("%s is not find" %file_name)
		return (myarray,myseq)
	mylines=[line.strip() for line in myfile.readlines()]
	myfile.close()
	for i in range (len(mylines)):
		if len(mylines[i]) >0 and mylines[i][0:4]=="Last":
			pssm_idx=i+2
			newscore=[]
			scoreline=mylines[pssm_idx].split()

			while len(scoreline)>0 and scoreline[0]!="Standard":
				myseq+=scoreline[1]
				for j in range(2,42):
					if j < 22:
						# for position-specific scoring matrix
						newscore.append(int(scoreline[j]))
					else:
						# for weighted observed percentages
						newscore.append(float(scoreline[j])/100)
				myarray.append(newscore)
				newscore=[]
				pssm_idx+=1
				scoreline=mylines[pssm_idx].split()
			break
	return (myarray,myseq)
#fr_n="../data/pssm.train/P31043.fasta.pssm"
#myarray,myseq=r_array(fr_n)
#print myarray
#print myseq
def conservation_features(query_array,seq):
	#compute the realtive entropy 
	AAdb={'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'X':20}
	#AAdb={'A':0,'D':1,'C':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19,'U':20,'X':20}

	r_ens=[[0 for j in range(20)] for i in range(20)] # add the scores

	t_ens=[[0 for j in range(20)] for i in range(20)] #normalized by length

	results=[[0 for j in range(20)] for i in range(20)]
	for i in range(len(query_array)):
		#relative_entropy=0.0
		eachline=query_array[i]
		try:
			idx=AAdb[seq[i]]
		except ValueError:
			idx=21
		if idx<20:
			for j in range(20):
				item=eachline[j]+r_ens[idx][j]
				r_ens[idx][j]=item
			#print "====="
			#print query_array[i]
			#print r_ens[idx]
			#r_ens[idx]=My_add_two_vectors(query_array[i],r_ens[idx])
			#r_ens[idx]=t_r_ens
	#normalized by the length 2011-04-18
	len_seq=len(seq)
	if len_seq > 0:
		for i in range(20):
			for j in range(20):
				t_ens[i][j]=r_ens[i][j]*1.0/len_seq
		#normalized by the length up 2011-04-18
		for i in range(20):
			#print r_ens[i]
			results[i]=My_pssm_normalize(t_ens[i])
			#print "+",results[i]

	return results
####
###---Main----
def get_idx(fr_n):
	""" read in fasta file; return idxs"""
	fr=open(fr_n,'r')
	lines=[line.strip() for line in fr.readlines()]
	fr.close()
	idxes=[]
	for i in range(len(lines)):
		if len(lines[i])>0 and lines[i][0]=='>':
			proid=lines[i][1:]
			idxes.append(proid)
	return idxes

def main(fx_n,folder,fw_n):

	# read in fasta file ; return protein id 
	fxlines=get_idx(fx_n)

	#fw_n="../data/lig.pssm.csv"#"../data/vol.pssm.csv"
	fw=open(fw_n,'w') #2010-08-15

	for i in range(len(fxlines)):
		if len(fxlines[i])>0 :
			proid=fxlines[i]
			f_name=folder+"/"+proid+".pssm" #2010-08-15
			myarray,fastaseq=r_array(f_name)
			if len(myarray)==0 or len(fastaseq)==0:
				print( "### %s is not find " %proid)
			con_features=[]
			con_features=conservation_features(myarray,fastaseq)
			#if len(myarray)==0 and len(fastaseq)==0:
			#	print con_features #check the protein, which has no pssm file
			for j in range(len(con_features)):
				for k in range(len(con_features[j])):
					fw.write("%.5f," %con_features[j][k])
			fw.write("\n")
	fw.close()

if __name__=="__main__":
	if len(sys.argv)<3:
		print("Usage: pssmfeature.py folderpath fr_n folder fw_n\nfr_n: fast; folder:where pssm; fw_n: pssm.csv")
	else:
		#folderpath=sys.argv[1]
		fx_n=sys.argv[1]
		folder=sys.argv[2]
		fw_n=sys.argv[3]
		main(fx_n,folder,fw_n)



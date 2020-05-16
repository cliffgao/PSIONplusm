#!/usr/bin/env python

""" compute the spinex"""

import os,sys
#sys.path.append("/home/cliff/gjz-cmd-line/")
from my_average import my_average
def mv_names():
	idx_f=open("d62.idx")
	idxs=[line.strip() for line in idx_f.readlines()]
	idx_f.close()

	for i in range(len(idxs)):
		fr_n=idxs[i]
		cmdline="mv ./profile/%s.pssm  ./profile/%s.mat" %(fr_n,fr_n)
		os.system(cmdline)
	return 1
#mv_names()
# compute the spinex
#cmdline2="../spX.pl  d62.idx profile/"
#os.system(cmdline2)

def ASA2RSA(AA,ASA,CUTOFF=0.25):
	"""1 for exposed; 0 for buried"""
	AXA={'A':110.2,'D':144.1,'C':140.4,'E': 174.7,'F': 200.7, \
'G':78.7, 'H':181.9,'I': 185.0 , 'K': 205.7, \
'L': 183.1,'M': 200.1, 'N':146.4,'P':141.9,'Q':178.6, \
'R':229.0,'S':117.2,'T':138.7 ,'V':153.7,\
'W':240.5,'Y':213.7}
	RSA=[]
	predicts=""
	if len(AA)!= len(ASA):
		print("length is different: AA:%d ASA:%d " %(len(AA),len(ASA)))
	else:
		for i in range(len(AA)):
			oneAsa=ASA[i]*1.0/AXA[AA[i]]
			RSA.append(oneAsa)
			if oneAsa > CUTOFF:
				onePred="1"
			else:
				onePred="0"
			predicts=predicts+onePred
	return (RSA,predicts)

def get_asa_spinex(fr_n):
	fr=open(fr_n,'r')
	predicts=[line.strip() for line in fr.readlines()]
	fr.close()
	ASA=[]
	AA=""
	SS=""
	for i in range(1, len(predicts)): # donot do with header
		if predicts[i]:
			onepredict=predicts[i].split()
			try:
				ASA.append(float(onepredict[10]))
				AA=AA+onepredict[1]
				SS=SS+onepredict[2]
			except IndexError:
				print("******Error: %s" %fr_n)
	return (AA,SS,ASA)




def my_dowith_spinex(fr_n,folder,fw_n,CUTOFF=0.25):
	"""input fasta,folder,fw_n, cutoff; ouput >id;seq;ss;exposed/burid"""
	fr=open(fr_n,'r')
	lines=[line.strip() for line in fr.readlines()]
	fr.close()
	AA=""
	SS=""
	ASA=RSA=[]

	#fw_n="d62.rsa.dat"   RSA values
	fw=open(fw_n,'w')

	for i in range(len(lines)):
		if len(lines[i])>0 and lines[i][0]=='>':
			proid=lines[i].split()[0][1:]
			proseq=lines[i+1]
			output_f=folder+ "/"+proid+".spXout" #"./spXout/%s.spXout" %proid.upper()
			AA,SS,ASA=get_asa_spinex(output_f)
			#print(ASA)
			if proseq!=AA:
				print("predicted vs. proseq different: %s" %proid)
				# do with the sequence contain X
				#write the id;proseq;
				fw.write("%s\n%s\n%s\n" %(lines[i],lines[i+1],SS))
				tmp_seq=proseq.replace("X","")
				tmp_RSA,tmp_BdEd=ASA2RSA(tmp_seq,ASA)
				#used the average rsa to replace the position x.
				tmp_ave_RSA=my_average(tmp_RSA)
				#write the buried/exposed
				count_x=0
				for k in range(len(proseq)):
					if proseq[k]=="X":
						count_x=count_x+1
						if tmp_ave_RSA> CUTOFF:
							fw.write("1")
						else:
							fw.write("0")
					else:
						if tmp_RSA[k-count_x] > CUTOFF:
							fw.write("1")
						else:
							fw.write("0")
				fw.write("\n")
				# write the rsa values
				count_x=0
				for k in range(len(proseq)):
					if proseq[k]=="X":
						fw.write("%.5f," %tmp_ave_RSA)
						count_x=count_x+1
					else:
						fw.write("%.5f," %tmp_RSA[k-count_x])
				fw.write("\n")
			else:
				RSA,BdEd=ASA2RSA(AA,ASA)
				#write the id; proseq;ss,rsa_binary
				fw.write("%s\n%s\n%s\n%s\n" %(lines[i],AA,SS,BdEd))
				# write the rsa values
				for j in range(len(RSA)):
					fw.write("%.5f," %RSA[j])
				fw.write("\n")
	fw.close()

if __name__=="__main__":
	if len(sys.argv)<3:
		print("Usage: my_spinex.py  fr_n  folder_name fw_n")
	else:
		fr_n=sys.argv[1]
		folder=sys.argv[2]
		fw_n=sys.argv[3]
		my_dowith_spinex(fr_n,folder,fw_n)
		#print("### Note RSA is compute based one cutoff 0.25\n>id;seq;ss;BdEd;rsa__value")

#./my_spinex.py lipp.represent.fa   rsa.dat  ./asa-all				


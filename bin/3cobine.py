#!/usr/bin/env python
"""
this script is used to combine the prediction of Spine3.0
AXA:
110.2 (Ala), 144.1 (Asp),140.4 (Cys), 174.7 (Glu), 200.7
(Phe), 78.7 (Gly), 181.9 (His), 185.0 (Ile), 205.7 (Lys),
183.1 (Leu), 200.1 (Met),146.4 (Asn), 141.9 (Pro), 178.6
(Gln), 229.0 (Arg), 117.2 (Ser), 138.7 (Thr), 153.7 (Val),
240.5 (Trp), and 213.7 (Tyr) respectively.
###just for CBTOPE  2010-11-23
>id
fastasequence
secondary structure
asa value
The sequence fragment 's 12w so large,
I need to search for the file
"""
import os,sys
#from get_idx import get_idx
def get_idx(fr_n):
	fr=open(fr_n,'r')
	lines=[line.strip() for line in fr.readlines()]
	fr.close()
	idxes=[]
	for i in range(len(lines)):
		if len(lines[i])>0 and lines[i][0]=='>':
			proid=lines[i][1:]
			idxes.append(proid)
	return idxes
def Fetch_file(name,flag):
  #name=P11387:22
  if flag=="ss":
    folder="./data/SS/"#folder_path #"../data/non.ssasa.results/"#"../data/lig.ssasa.results/"
  elif flag=="asa":
    folder="./data/ASA/"#folder_path #"../data/non.ssasa.results/"
  for i in range(1):
    filename1=folder+name+"."+flag
    lines=[]
    #filename2=folder+str(i+1)+"/"+name.replace(":","_")+"."+flag
    if os.path.isfile(filename1):
      fr=open(filename1,'r')
      lines=[line.strip() for line in fr.readlines()]
      fr.close()
      break
    #elif os.path.isfile(filename2):
    #  fr2=open(filename2,'r')
    #  lines=[line.strip() for line in fr2.readlines()]
    #  fr2.close()
    #  break
    else:
      lines=None
  return lines
#lines=Fetch_file("P14917:1","asa")
#print lines[0:4]

#"""
def main(idx_f,fw_n):

	AXA={'A':110.2,'D':144.1,'C':140.4,'E': 174.7,'F': 200.7, \
'G':78.7, 'H':181.9,'I': 185.0 , 'K': 205.7, \
'L': 183.1,'M': 200.1, 'N':146.4,'P':141.9,'Q':178.6, \
'R':229.0,'S':117.2,'T':138.7 ,'V':153.7,\
'W':240.5,'Y':213.7}
	
	#fr=open(idx_f,'r')#"../data/epitopia.idx",'r')
#fr=open("../structure_data/substr20.struct.idx",'r')
#fr=open("../structure_data/new_fragment_3D.idx",'r')
	#idxs=[line.strip() for line in fr.readlines()]
	#fr.close()
	idxs=get_idx(idx_f)
	f3=open(fw_n,'w') #"../data/non.ss.rsa.txt",'w')
#f3=open("../structure_data/epitopia.stru.3d.ss.rsa.txt",'w')
#folder_name="../structure_data/struct_3D_fragment_pssm_ssasa"#"../structure_data/struct-fragment-ssasa" #"../epitopia.ssasa.results"

	for i in range(len(idxs)): # modify this 
  # read in the secdary structure
  #proid="chen.p.%d" %(i)
  #proid=idxs[i].replace(":","_")
		proid=idxs[i]
		f3.write(">%s\n" %proid)

  # secondary strcuture
  #ss_file_name="%s/%s.ss" %(folder_name,idxs[i])
  #if not os.path.isfile(ss_file_name): #---08-14
  #  print "%s is not exist" %idxs[i]
  #else:
  #  ss_file     = open(ss_file_name,'r')
  #  sslines     = [line.strip() for line in ss_file.readlines()]
  #  ss_file.close()  #----o8-14
		sslines=Fetch_file(proid,"ss")
		asalines=Fetch_file(proid,"asa")
		if (sslines is None) or (asalines is None):
			print( "%s is not exist" %proid)
		else:
  #get ss sequence
			for j in range(len(sslines)):
				if len(sslines[j])>0 :
					if sslines[j]=="[sequence]":
						f3.write("%s\n" %sslines[j+1])
					if sslines[j]=="[secondary]":
						f3.write("%s\n" %sslines[j+1])
  # read in the solvent structure
    #asa_file_name="%s/%s.asa" %(folder_name,idxs[i])
    #asa_file    = open(asa_file_name,'r')
    #asalines   = [line.strip() for line in asa_file.readlines()]
    #asa_file.close()
			for k in range (2,len(asalines)):
				if len(asalines[k])>0:
					asaeach=asalines[k].split()
      #write the RSA value
					f3.write("%.4f," %(float(asaeach[2])/AXA[asaeach[1]]))
			f3.write("\n")

	f3.close()    
#  print ss_file_name
#"""
if __name__=="__main__":
	if len(sys.argv)<2:
		print("Usage: 1cobine.py idx_f fw_n")
	else:
		idx_f=sys.argv[1]
		fw_n=sys.argv[2]
		#folder_path="./ASASS/"#sys.argv[3]

		main(idx_f, fw_n)
 

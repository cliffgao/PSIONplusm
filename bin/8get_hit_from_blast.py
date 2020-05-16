#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
import sys

class Blast:
	def __init__(self,fr_n):
		self.fr_n=fr_n
	def get_content(self):
		fr=open(self.fr_n,'r')
		lines=[line.strip() for line in fr.readlines()]
		fr.close()
		return lines

	def get_result(self,fw_n,eCutoff=10):
		lines=self.get_content()
		fw=open(fw_n,'w')
		for i in range(len(lines)):
			if lines[i]:
				eachline=lines[i].split()
				if eachline[0]=="Query=":
					fw.write(">%s\n" %(eachline[1]))
					#print (eachline[1])
				#if eachline[0]=="Searching.done":
				if lines[i].split(".")[0]=="Searching":
					#print ("OK")
					k=i+2
					
					nextline=lines[k]
					while(len(nextline)==0):
						k=k+1
						nextline=lines[k]
					while(nextline):
						if nextline=="***** No hits found ******":
							fw.write("--\n")
							break
						if nextline=="Sequences producing significant alignments:                      (bits) Value":
							#hit_ids=lines[k+2].split()
							#hit_id=hit_ids[0]
							#e_value=float(hit_ids[-1])
							#print("%.3f" %(float(hit_ids[-1])))
							#fw.write("%s\n" %hit_id)
							u=k+2
							hitflag=False
							while(len(lines[u])>0 and lines[u][0]!='>'):
								if len(lines[u])>0:
									hit_ids=lines[u].split()
									hit_id=hit_ids[0]
									e_valuex=hit_ids[-1]
									if e_valuex[0]=='e':
										e_value="1"+e_valuex
									else:
										e_value=e_valuex
									#print(hit_id,e_value)
									if float(e_value) < eCutoff:
										hitflag=True
										fw.write("%s  %s\n" %(hit_id,e_value))
										
								u=u+1
							k=u
						#if nextline[0]=='>':
							#break
							if hitflag==False:
								fw.write("--\n")
							break

						k=k+1
						nextline=lines[k]       
		fw.close()

def main(fr_n,fw_n,eCutoff):
	#fr_n="./blast.Types.out"
	myblast=Blast(fr_n)
	#fw_n="./eg.types.fa"
	myblast.get_result(fw_n,eCutoff)

	#fr_n="./myblast.e1.out"
	#myblast=Blast(fr_n)
	#fw_n="./hit.e1.txt"
	#myblast.get_result(fw_n)
#fr_n="./blast.Types.out"
#myblast=Blast(fr_n)
#fw_n="./eg.types.fa"
#myblast.get_result(fw_n,0.005)


#"""
if __name__=="__main__":
	if len(sys.argv)<3:
		print("Usage: blast.results.py  blast.out   myblast.fa")
	else:
		fr_n=sys.argv[1]
		fw_n=sys.argv[2]
		eCutoff=float(sys.argv[3])
		main(fr_n,fw_n,eCutoff)
#"""

#!/usr/bin/env python

import sys
def csv2svm(fr_n):
	""" change the csv file into svm format
	label is the first one"""
	fr=open(fr_n,'r')
	lines=[line.strip() for line in fr.readlines()]
	fr.close()

	fw_n=fr_n+".svm"

	fw=open(fw_n,'w')

	for i in range(len(lines)):
		if len(lines[i])>0 :
			eachline=lines[i].split(",")

			fw.write("%s " %eachline[0]) # write the head
			svm_idx=1
			for j in range(1,len(eachline)-1): #rm the last ","
				fw.write("%d:%s " %(svm_idx,eachline[j]))
				svm_idx=svm_idx+1
			fw.write("\n")
	fw.close()

	return 1

if __name__=="__main__":
	if len(sys.argv)<2:
		print("Usage: csv2svm.py  fr_n \n The fr_n.svm is the results file; Head is the label")
	else:
		fr_n=sys.argv[1]
		csv2svm(fr_n)


#!/usr/bin/env python
#""" get the not-actual  label for test dataset """
import sys

def get_label(fr_n,fw_n):
	fr=open(fr_n,'r')
	lines=[line.strip() for line in fr.readlines()]
	fr.close()

	fw=open(fw_n,'w')
	for i in range(len(lines)):
		if len(lines[i]) > 0 and lines[i][0]=='>':
			fw.write("-1\n") #print tmp label 
	fw.close()
	return 1	

fr_n=sys.argv[1]
fw_n=sys.argv[2]
#flag=sys.argv[3]
get_label(fr_n,fw_n)

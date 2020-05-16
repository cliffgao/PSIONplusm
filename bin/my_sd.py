#!/usr/bin/env python

""" compute the standard deviation"""
import sys

def my_sd(mylistx):
	"""compute the standard deviation"""
	mylist=[float(item) for item in mylistx]
	N=len(mylist)
	if N>0:
		average=sum(mylist)/N
		sd=0.0
		for i in range(N):
			sd+= (mylist[i]-average)**2
		return (sd/N)**0.5
	else:
		print("length of list is zero")
		return 0
def my_sd2(mylistx):
	"""compute the standard deviation"""
	mylist=[float(item) for item in mylistx]
	N=len(mylist)
	if N>0:
		average=sum(mylist)/N
		sd=0.0
		for i in range(N):
			sd+= (mylist[i]-average)**2
		return (sd/(N-1))**0.5
	else:
		print("length of list is zero")
		return 0
if __name__=="__main__":
	if len(sys.argv)<2:
		print("usage: my_sd.py mylist")
	else:
		my_sd(mylist)

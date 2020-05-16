#!/usr/bin/env python
"""
return the average/mean of mylist
"""


def my_average(mylist):
	size=len(mylist)
	results=0.0
	#print "********",mylist
	if None in mylist:
		print(" mylist contains None")
		print( mylist)
		return results
	elif size==0:
		results=0.0
	else:
		results=sum(mylist)*1.0/size

	return results

#
#print my_average([1,2,3,4])
      

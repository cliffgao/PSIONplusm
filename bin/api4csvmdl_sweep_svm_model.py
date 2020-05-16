#!/usr/bin/env python 

"""
select the features according the select feature index
this script is modify in 2010-08-05 
for just modifying just one file

"""
import math
import os
import sys
def Get_value(mylist,idx):
    value=0
    for i in range(len(mylist)):
        item=int(mylist[i].split(":")[0])
        if idx==item:
            value=mylist[i].split(":")[1]
            break
    
    return value

def Get_features_order(fr_n,flagType):
    fr=open(fr_n,'r')
    lines=[line.strip() for line in fr.readlines()]
    fr.close()
    results=[]
    if flagType=="row":
        for i in range(len(lines)):
            if len(lines[i])>0 and lines[i][0:6]=="select":
                t_m=lines[i].split("[")[1].split("]")[0]
                t_m2=t_m.split(",")
                for j in range(len(t_m2)):
                    results.append(int(t_m2[j]))
    elif flagType=="col":
        for i in range(len(lines)):
            if len(lines[i])>0 :
                results.append(int(lines[i]))
    #results=[int(line.strip()) for line in lines[0].split(",")]
    # return the features idx
    return results
#print Get_features_order("vol.lig.model.svm.select")
#"""

def Modify_file(infile,outfile,features_n,flagType):
    ##infile, original file; outfile, output file; 
    #feature_order, the feature 
    #
    fr=open(infile,'r')
    ori_contant=[line.strip() for line in fr.readlines()]
    fr.close()
    # write file
    fw=open(outfile,'w')
    #get features order
    #
    feature_order=Get_features_order(features_n,flagType)
    #print len(feature_order),feature_order

    for j in range(len(ori_contant)):
        each_contant=ori_contant[j].split()
        fw.write("%s " %each_contant[0]) # write the label
        svm_idx=1
        for k in range(len(feature_order)):
            idx=int(feature_order[k])
            search_content=each_contant[1:]
            value=Get_value(search_content,idx)
            if value=="noValue":
                print(" idx: %d  is out of range %d " %(idx,len(each_contant)))
            else:
                fw.write("%d:%s " %(svm_idx,value) )
                svm_idx+=1
        fw.write("\n")
    # write the test file
    fw.close()
    
    return 0 

#------------------

if __name__=="__main__":
    if len(sys.argv)<4:
        print("Usage: sweep_svm_model.py original_svm_file output_file  features_idx_file flagType(row,col?)")
    else:
        svm_model_n=sys.argv[1]
        out_model_n=sys.argv[2]
        features_idx_n=sys.argv[3]
        flagType=sys.argv[4]
        Modify_file(svm_model_n,out_model_n,features_idx_n,flagType)
#"""

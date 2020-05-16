#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 18:47:52 2019

@author: Cliff
"""
from csvmdl_class import CSVMDL
from PSION_features_class import PSION_features
import os,sys
import numpy as np 
def output_svm_mdl(frn,scale_info_fn,selectedFn,svmdl):
    
    fr=open(frn,'r')
    lines=[line.strip() for line in fr.readlines()]
    fr.close()
    
    proseq=lines[1] #sys.argv[1]
    proid=lines[0].split()[0][1:] #sys.argv[2]
    ### get the label and all features
    cls_fs=PSION_features(proseq,dir)
    ##labels and all features
    xfull_lb_fs=cls_fs.get_label_features()
    ### change .5f the same as the original model
    full_lb_fs=[]
    full_lb_fs.append(xfull_lb_fs[0])
    for k in range(1,len(xfull_lb_fs)):
        item=round(xfull_lb_fs[k],5)
        full_lb_fs.append(item)
    ### scaled features and selected optimized features
    cls_scale_fs=CSVMDL(full_lb_fs,scale_info_fn)
    cls_scale_fs.get_svm(selectedFn,svmdl)
    
    #print("svm model is %s" %svmdl)
    ### run svm
    return(0)

if __name__=="__main__":
    if len(sys.argv)<5:
        print("Usage fastafn, scale_info_fn, selectedFn,svmdl")
        print(sys.argv[1])
        print(sys.argv[2])
        print(sys.argv[3])
        print(sys.argv[4])
    else:
        frn=sys.argv[1]
        dir=sys.argv[2]
        scale_info_fn=sys.argv[3]
        selectedFn=sys.argv[4]
        svmdl=sys.argv[5]
        output_svm_mdl(frn,scale_info_fn,selectedFn,svmdl)

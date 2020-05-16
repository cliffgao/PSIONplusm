#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 19:05:22 2019

@author: Cliff
"""
from api4csvmdl_scale_svmodel_test import get_scale_file
from api4csvmdl_sweep_svm_model import Get_features_order
class CSVMDL:
    def __init__(self,LabelFeatures,scaleInfoFn):
        ###
        self.mymin=[]
        self.mymax=[]
        #scaleInfoFn="./data/
        mymin,mymax=get_scale_file(scaleInfoFn)
        self.mymin=mymin
        self.mymax=mymax
        self.LabelFeatures=LabelFeatures
        
    def scale_model_according_minmax(self):
        ### scaled features accroding to minmax
        #EPS=0.00000001  
        features=self.LabelFeatures[1:] 
        label=self.LabelFeatures[0]
        item_ma=self.mymax
        item_mi=self.mymin
        scaled_fs=[]
        scaled_lb_fs=[]
        #    #scale the features
        for j in range(len(features)):
            if item_ma[j]==item_mi[j] :  # if the fenmu is zero
                item=item_ma[j]
            else:
                item=(features[j]*2.0-(item_mi[j]+item_ma[j]))/(item_ma[j]-item_mi[j])
            scaled_fs.append(item)
# =============================================================================
#             myrange=item_ma[j]-item_mi[j]
#             if abs(myrange)<EPS:  # if the fenmu is zero
#                 item=item_ma[j]
#             else:
#                 item=(features[j]*2.0-(item_mi[j]+item_ma[j]))/(myrange)
#             #print("%d:%.5f " %(svm_idx,item)),
# =============================================================================
        scaled_lb_fs=[label]+scaled_fs
        
        return (scaled_lb_fs)
    def get_scaled_lb_fs(self):
        return(self.scale_model_according_minmax())

    def sweep_features_according_fn(self,selectedFn):
        flagType="col"
        #begin with 1
        selectedIdx=Get_features_order(selectedFn,flagType)
        selected_lb_fs=[]
        
        scaled_lb_fs=self.scale_model_according_minmax()
        scaled_lb=scaled_lb_fs[0]

        selected_lb_fs.append(scaled_lb)
        for aidx in selectedIdx:
            item=scaled_lb_fs[aidx] #aidx begin with 1; 
            #scaled_lb_fs  begin with label
            selected_lb_fs.append(item)
        return(selected_lb_fs)
        
    def csv2svm(self,labelFeatures,fwn):
    #change the csv file into svm format
    #label is the first one"""

        fw=open(fwn,'w')
    
        eachline=labelFeatures #1st is label
        fw.write("%s " %eachline[0]) # write the head
        svm_idx=1
        for j in range(1,len(eachline)): #rm the last ","
            fw.write("%d:%.4f " %(svm_idx,eachline[j]))
            svm_idx=svm_idx+1
        fw.write("\n")
        fw.close()
        return(fwn)
    def get_svm(self,selectedFn,fwn):
        selected_lb_fs=self.sweep_features_according_fn(selectedFn)
        self.csv2svm(selected_lb_fs,fwn)
        

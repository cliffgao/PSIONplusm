#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 10:55:41 2019

@author: Cliff
nk_pi@hotmail.com
"""
import os,sys 
from api_01_AA49 import api_aacv_pp49
from api_02_dowith_psipred import api_dowith_psipred
from api_02_ss_features import api_ss_features
from api_03_dowith_disopred import api_dowith_disorder
from api_03_dis_features import api_dis_features
from api_04_dowith_spinex import api_dowith_spinex
from api_04_rsa_features import api_rsa_features
from api_05_pssm_feature import  api_pssm_features
from api_06_diAA import api_my_diAA

class PSION:
    ### PSION is a model build by SVM
    def __init__(self,proseq,proid="QuerySeq"):
# =============================================================================
#         # read fasta file
#         fr=open(frn,'r')
#         lines=[line.strip() for line in fr.readlines()]
#         fr.close()
#         ###
# =============================================================================
        self.faacvpp49=[]
        self.fss=[]
        self.fdis=[]
        self.frsa=[]
        self.fcons=[]
        self.fdiaa=[]
        psipredfn="./data/ss.test/%s.ss" %proid
        disorderfn="./data/dis.test/%s.diso"%proid
        spinexfn="./data/spx.test/%s.spXout" %proid
        pssmfn="./data/pssm.test/%s.pssm" %proid
# =============================================================================
#         ### for each fasta id
#         for i in range(len(lines)):
#             if len(lines[i])>0 and lines[i][0]=='>':
#                 proid=lines[i].split()[0][1:]
#                 proseq=lines[i+1]
       
        ### aa composition+ 49 properties
        self.faacvpp49=api_aacv_pp49(proseq)
        ### ss from psipred
        if os.path.isfile(psipredfn):
            SS=api_dowith_psipred(psipredfn)
        else:
            SS=""
        self.fss=api_ss_features(proseq,SS)
        ### disorder prediction
        if os.path.isfile(disorderfn):
            DIS=api_dowith_disorder(disorderfn)
        else:
            DIS=""
        self.fdis=api_dis_features(proseq,DIS)
        ### rsa features
        if os.path.isfile(spinexfn):
            spdAA,spdSS,spdBdEd,spdRSA=api_dowith_spinex(spinexfn,proseq)
        else:
            spdAA=""
            spdSS=""
            spdBdEd=""
            spdRSA=[]
        self.frsa=api_rsa_features(spdAA,spdSS,spdRSA)
        ### ### conservation scores
        if os.path.isfile(pssmfn):
            self.fcons=api_pssm_features(pssmfn,proseq)
        else:
            self.fcons=[]
        ### diaa    
        self.fdiaa=api_my_diAA(proseq)
    def get_aapp49_features(self):
        return(self.faacvpp49)
    def get_ss_features(self):
        return(self.fss)
    def get_dis_features(self):
        return(self.fdis)
    def get_rsa_features(self):
        return(self.frsa)
    def get_con_features(self):
        return(self.fcons)
    def get_diaa_features(self):
        return(self.fdiaa)
    def get_all_features(self):
        all_features=[]
        for item in self.faacvpp49:
            all_features.append(item)
        for item in self.fss:
            all_features.append(item)
        for item in self.fdis:
            all_features.append(item)
        for item in self.frsa:
            all_features.append(item)
        for item in self.fcons:
            all_features.append(item)
        for item in self.fdiaa:
            all_features.append(item)
        return(all_features)
    def get_label_features(self):
        pesudoLabel=[-1]
        all_fs=self.get_all_features()
        mdl=pesudoLabel+all_fs
        ## labels+allfeatures
        return(mdl)

if __name__=="__main__":
    if len(sys.argv)<2:
        print("Usage proseq proid")
    else:
        frn=sys.argv[1]
        fr=open(frn,'r')
        lines=[line.strip() for line in fr.readlines()]
        fr.close()
        
        proseq=lines[1] #sys.argv[1]
        proid=lines[0].split()[0][1:] #sys.argv[2]
        aa=PSION(proseq,proid)
        print(aa.get_label_features())
        

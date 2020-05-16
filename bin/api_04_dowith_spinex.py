#!/usr/bin/env python

""" compute the spinex"""

import os,sys
#sys.path.append("/home/cliff/gjz-cmd-line/")
from my_average import my_average
def mv_names():
    idx_f=open("d62.idx")
    idxs=[line.strip() for line in idx_f.readlines()]
    idx_f.close()

    for i in range(len(idxs)):
        fr_n=idxs[i]
        cmdline="mv ./profile/%s.pssm  ./profile/%s.mat" %(fr_n,fr_n)
        os.system(cmdline)
    return 1
#mv_names()
# compute the spinex
#cmdline2="../spX.pl  d62.idx profile/"
#os.system(cmdline2)

def ASA2RSA(AA,ASA,CUTOFF=0.25):
    """1 for exposed; 0 for buried"""
    AXA={'A':110.2,'D':144.1,'C':140.4,'E': 174.7,'F': 200.7, \
'G':78.7, 'H':181.9,'I': 185.0 , 'K': 205.7, \
'L': 183.1,'M': 200.1, 'N':146.4,'P':141.9,'Q':178.6, \
'R':229.0,'S':117.2,'T':138.7 ,'V':153.7,\
'W':240.5,'Y':213.7}
    RSA=[]
    predicts=""
    if len(AA)!= len(ASA):
        print("length is different: AA:%d ASA:%d " %(len(AA),len(ASA)))
    else:
        for i in range(len(AA)):
            oneAsa=ASA[i]*1.0/AXA[AA[i]]
            RSA.append(oneAsa)
            if oneAsa > CUTOFF:
                onePred="1"
            else:
                onePred="0"
            predicts=predicts+onePred
    return (RSA,predicts)

def get_asa_spinex(fr_n):
    fr=open(fr_n,'r')
    predicts=[line.strip() for line in fr.readlines()]
    fr.close()
    ASA=[]
    AA=""
    SS=""
    for i in range(1, len(predicts)): # donot do with header
        if predicts[i]:
            onepredict=predicts[i].split()
            try:
                ASA.append(float(onepredict[10]))
                AA=AA+onepredict[1]
                SS=SS+onepredict[2]
            except IndexError:
                print("******Error: %s" %fr_n)
    return (AA,SS,ASA)




def api_dowith_spinex(output_f,proseq,CUTOFF=0.25):
    """input fasta,folder,cutoff; ouput >id;seq;ss;exposed/burid"""
    AA,SS="",""
    ASA=[]
    RSA=[]
    BdEd=""
    if os.path.isfile(output_f):
        AA,SS,ASA=get_asa_spinex(output_f)
    else:
        print("Can not find SpineX prediction of %s" %output_f)
        
    #print(ASA)
    if proseq!=AA:
        print("predicted vs. proseq different: %s" %output_f)
        # do with the sequence contain X
        #write the id;proseq;
        tmp_seq=proseq.replace("X","")
        tmp_RSA,tmp_BdEd=ASA2RSA(tmp_seq,ASA)
        #used the average rsa to replace the position x.
        tmp_ave_RSA=my_average(tmp_RSA)
        #write the buried/exposed
        count_x=0
        for k in range(len(proseq)):
            if proseq[k]=="X":
                count_x=count_x+1
                if tmp_ave_RSA> CUTOFF:
                    BdEd+="1" #fw.write("1")
                else:
                    BdEd+="0" #fw.write("0")
            else:
                if tmp_RSA[k-count_x] > CUTOFF:
                    BdEd+="1"#fw.write("1")
                else:
                    BdEd+="0" #fw.write("0")
        print("Meet %d X" %count_x)
        # write the rsa values
        for k in range(len(proseq)):
            if proseq[k]=="X":
                RSA.append(tmp_ave_RSA)
            else:
                RSA.append(tmp_RSA[k-count_x])
    else:
        RSA,BdEd=ASA2RSA(AA,ASA)
           
    return((AA,SS,BdEd,RSA))


def my_dowith_spinex(fr_n,folder,CUTOFF=0.25):
    """input fasta,folder,cutoff; ouput >id;seq;ss;exposed/burid"""
    fr=open(fr_n,'r')
    lines=[line.strip() for line in fr.readlines()]
    fr.close()
   
    SSS=[]
    BdEdS=[]
    RSAS=[]

    for i in range(len(lines)):
        if len(lines[i])>0 and lines[i][0]=='>':
            proid=lines[i].split()[0][1:]
            proseq=lines[i+1]
            BdEd=""
            RSA=[]
            #output_f=folder+ "/"+proid+".spXout" #"./spXout/%s.spXout" %proid.upper()
            output_f="%s/%s.spXout" %(folder,proid)
            if os.path.isfile(output_f):
                AA,SS,ASA=get_asa_spinex(output_f)
            else:
                AA,SS="",""
                ASA=[]
                print("Can not find SpineX prediction of %s" %proid)
                
            #print(ASA)
            if proseq!=AA:
                print("predicted vs. proseq different: %s" %proid)
                # do with the sequence contain X
                #write the id;proseq;
                tmp_seq=proseq.replace("X","")
                tmp_RSA,tmp_BdEd=ASA2RSA(tmp_seq,ASA)
                #used the average rsa to replace the position x.
                tmp_ave_RSA=my_average(tmp_RSA)
                #write the buried/exposed
                count_x=0
                for k in range(len(proseq)):
                    if proseq[k]=="X":
                        count_x=count_x+1
                        if tmp_ave_RSA> CUTOFF:
                            BdEd+="1" #fw.write("1")
                        else:
                            BdEd+="0" #fw.write("0")
                    else:
                        if tmp_RSA[k-count_x] > CUTOFF:
                            BdEd+="1"#fw.write("1")
                        else:
                            BdEd+="0" #fw.write("0")
                print("Meet %d X" %count_x)
                # write the rsa values
                for k in range(len(proseq)):
                    if proseq[k]=="X":
                        RSA.append(tmp_ave_RSA)
                    else:
                        RSA.append(tmp_RSA[k-count_x])
            else:
                RSA,BdEd=ASA2RSA(AA,ASA)
            SSS.append(SS)
            BdEdS.append(BdEd)
            RSAS.append(RSA)
            
    return((AA,SS,BdEdS,RSAS))    
    
    
# =============================================================================
# if __name__=="__main__":
#     if len(sys.argv)<3:
#         print("Usage: my_spinex.py  fr_n  folder_name fw_n")
#     else:
#         fr_n=sys.argv[1]
#         folder=sys.argv[2]
#         AA,SS,BdEdS,RSAS=my_dowith_spinex(fr_n,folder)
#         print(AA)
#         print(SS)
#         print(BdEdS)
#         print(RSAS)
# =============================================================================
 

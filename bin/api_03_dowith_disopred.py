#!/usr/bin/env python

""" compute the spinex"""

import os,sys

def do_with_disopred(fr_n):
    fr=open(fr_n,'r')
    predicts=[line.strip() for line in fr.readlines()]
    fr.close()
    
    AA=""
    DIS=""
    start=0
    #find the start
    for i in range(len(predicts)):
        if len(predicts[i])>0 and predicts[i].split()[0]=="1":
            start=i
            
    for i in range(start,len(predicts)): # donot do with header
        if predicts[i]:
            onepredict=predicts[i].split()
            DIS=DIS+onepredict[2]
            AA=AA+onepredict[1]
    #* disorder,.order
    return (AA,DIS)




def my_dowith_disopred(fasta_f,folder):
    """read in fasta;folder;fw_n"""
    fr_n=fasta_f
    fr=open(fr_n,'r')
    lines=[line.strip() for line in fr.readlines()]
    fr.close()
    AA=""
    DIS=""
    AAs=[]
    DISs=[]
    for i in range(len(lines)):
        if len(lines[i])>0 and lines[i][0]=='>':
            proid=lines[i].split()[0][1:]
            #proseq=lines[i+1]
            
            output_f="%s/%s.diso" %(folder,proid)
            if os.path.isfile(output_f):
                AA,DIS=do_with_disopred(output_f)
            else:
                AA=""
                DIS=""
                print("Not find DisoPred prediction of %s" %proid)
            AAs.append(AA)
            DISs.append(DIS)
    return((AAs,DISs))

def api_dowith_disorder(output_f):
    ## return predicted disorder
    AA=""
    DIS=""
    if os.path.isfile(output_f):
        AA,DIS=do_with_disopred(output_f)
    else:
        AA=""
        DIS=""
        print("Not find DisoPred prediction of %s" %output_f)
    return(DIS)
    
# =============================================================================
# if __name__=="__main__":
#     if len(sys.argv)<3:
#         print("Usage: my_dowith_psipred.py fasta_f, folder, fw_n")
#     else:
#         fasta_f=sys.argv[1]
#         folder=sys.argv[2]
#         
#         aas,diss=my_dowith_disopred(fasta_f,folder)
#         print(aas)
#         print(diss)
#                 
# 
# =============================================================================

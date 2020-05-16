#!/usr/bin/env python

""" compute the spinex"""

import os,sys

def do_with_psipred(fr_n):
    fr=open(fr_n,'r')
    predicts=[line.strip() for line in fr.readlines()]
    fr.close()
    
    AA=""
    SS=""
    for i in range(len(predicts)): # donot do with header
        if predicts[i]:
            onepredict=predicts[i].split()
            SS=SS+onepredict[2]
            AA=AA+onepredict[1]
    return (AA,SS)




def my_dowith_psipred(fasta_f,ssfolder):
    """read in fasta;psipred output folder"""
    fr_n=fasta_f
    fr=open(fr_n,'r')
    lines=[line.strip() for line in fr.readlines()]
    fr.close()
    AA=""
    SS=""
    AAS=[]
    SSS=[]
    for i in range(len(lines)):
        if len(lines[i])>0 and lines[i][0]=='>':
            proid=lines[i].split()[0][1:]
            #proseq=lines[i+1]
            output_f="%s/%s.ss" %(ssfolder,proid)
            if os.path.isfile(output_f):
                AA,SS=do_with_psipred(output_f)
            else:
                AA=""
                SS=""
                print("can not find PSIPRED prediction files %s" %proid)
            AAS.append(AA)
            SSS.append(SS)
    return((AAS,SSS))
    
def api_dowith_psipred(psipredfn):
    ### return predicted ss
    AA=""
    SS=""
    if os.path.isfile(psipredfn):
        AA,SS=do_with_psipred(psipredfn)
    else:
        AA=""
        SS=""
        print("can not find PSIPRED prediction files %s" %psipredfn)
    return(SS)
       
# =============================================================================
# if __name__=="__main__":
#     if len(sys.argv)<3:
#         print("Usage: my_dowith_psipred.py fasta_f, folder, fw_n")
#     else:
#         fasta_f=sys.argv[1]
#         folder=sys.argv[2]
#         aas,sss=my_dowith_psipred(fasta_f,folder)
#         print(aas)
#         print(sss)
#                 
# =============================================================================


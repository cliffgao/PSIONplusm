#!/usr/bin/env python
'''
this script is used to compute the
di-peptide 20*20=400
'''
import sys

def api_my_diAA(seq):
    ###
    AAdb={'A':0 ,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8, 'I':9,
            'L':10, 'K':11 ,'M':12,  'F':13,  'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'U':20,'X':20}
    ###
    seqlen=len(seq)
    aDIAA=[]
    #sigleAA=[0 for col in range(20)]
    #store the residue pair 
    doubleAA=[0 for col in range(400)]
    for j in range (seqlen-1):
        if ( AAdb[seq[j]]<20 and AAdb[seq[j+1]]<20):
            #idxI=AAdb[seq[j]]
            idxII=AAdb[seq[j]]*20 +AAdb[seq[j+1]]
            #sigleAA[idxI]=sigleAA[idxI]+1
            doubleAA[idxII]=doubleAA[idxII]+1
    #    e    ntropyIsum=0
    IIsum=sum(doubleAA)
    for k in range( 400):
        percentage=doubleAA[k]*1.0/IIsum
        aDIAA.append(percentage*100)
    #  each item 400
    return (aDIAA)


def diAA(fr_n):
    fastafile=open(fr_n,'r')
    fastalines=[line.strip() for line in fastafile.readlines()]
    fastafile.close()
    AAdb={'A':0 ,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8, 'I':9,
            'L':10, 'K':11 ,'M':12,  'F':13,  'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'U':20,'X':20}
    ###
    DIAAS=[]
    for i in range (len ( fastalines)):
        if (len(fastalines[i])>0 and fastalines[i][0]=='>'):
            seq=fastalines[i+1]
            seqlen=len(seq)
            aDIAA=[]
            #sigleAA=[0 for col in range(20)]
            #store the residue pair 
            doubleAA=[0 for col in range(400)]
            for j in range (seqlen-1):
                if ( AAdb[seq[j]]<20 and AAdb[seq[j+1]]<20):
                    #idxI=AAdb[seq[j]]
                    idxII=AAdb[seq[j]]*20 +AAdb[seq[j+1]]
                    #sigleAA[idxI]=sigleAA[idxI]+1
                    doubleAA[idxII]=doubleAA[idxII]+1
            #    e    ntropyIsum=0
            IIsum=sum(doubleAA)
            for k in range( 400):
                percentage=doubleAA[k]*1.0/IIsum
                aDIAA.append(percentage*100)
            DIAAS.append(aDIAA)
    # 2d  each item 400
    return (DIAAS)

# =============================================================================
# if __name__=="__main__":
#     if len(sys.argv)<2:
#         print("Usage diAA.py fr_n fw_n")
#     else:
#         fr_n=sys.argv[1]
#         #fw_n=sys.argv[2]
#         DIAAS=diAA(fr_n)
#         print(DIAAS)
# 
# =============================================================================


    

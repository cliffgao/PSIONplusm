#!/usr/bin/env python
"""
conbine the AAcontent and amino acid properity  # add the mean_var_49 in 2012-02-24
"""
import sys
#sys.path.append("/home/cliff/gjz-cmd-line/")
from my_sd import my_sd2

def cvAA(proseq):
    AAdb={'A':0 ,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10, 'K':11 ,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'U':20,'X':20}
    aa_num=[0 for i in range(20)] # count the number
    freq=[0 for i in range(20)] # frequence
    for i in range(len(proseq)):
        idx=AAdb[proseq[i]]
        if idx<20:
            aa_num[idx]+=1
    total=sum(aa_num)
    if total!=0:
        for j in range(20):
            freq[j]=aa_num[j]*1.0/total
    else:
        print("### total number is zero")
    return (freq)

def mean_var_49(proseq):
    ### amino acid property 
    dicfile=open("./bin/property.fa.normal",'r')
    diclines=[line.strip() for line in dicfile.readlines()]
    dicfile.close()

    results=[] #(mean_1,var_1,..., mean_49,var_49)
    # generate the dictionary
    AAindex=[]
    ss="ARNDCQEGHILKMFPSTWYV"
    for i in range(len(diclines)):
        if len(diclines[i])>0 and diclines[i][0]!='>':
            eachline=diclines[i].split()
            mydict={}.fromkeys(('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','U','X'))
            for j in range(len(eachline)):
                mydict[ss[j]]=float(eachline[j])
                mydict['U']=0.0
                mydict['X']=0.0
            AAindex.append(mydict)
    
    seq=proseq
    seqlen=len(seq)
    for j in range(len(AAindex)):# j-th dictionary
        scales=[]
        for k in range(seqlen):
            item_value=AAindex[j][seq[k]]
            scales.append(item_value)
        #compute average and stardard variance
        average=sum(scales)/(seqlen) # compute the average
        myvar=my_sd2(scales)                
        results.append(average)
        results.append(myvar)
    # compute average,standard deviation for each property.
    return (results)


def  api_aacv_pp49(proseq):
    # read in fasta sequence and return 20 amino aicd compostion; 
    # average, standard deviation of 49 properties.
    #compute the amino acid composition
    #AACV_mean_vars=[] #20 amino acid composition+mean/var of 49 properties
    
    aacontent=cvAA(proseq)
    #compute the k-property average,var
    proper49=mean_var_49(proseq)
    aacvpp=aacontent+proper49
    #AACV_mean_vars.append(aacontent+proper49)
    # one id one list (20+49*2=118) 
    # 2d list
    return(aacvpp)


def  api_aacv_pp49_fr(fr_n):
    # read in fasta sequence and return 20 amino aicd compostion; 
    # average, standard deviation of 49 properties.
    fr=open(fr_n,'r')
    lines=[line.strip() for line in fr.readlines()]
    fr.close()
    #fw=open(fw_n,'w')
    AACV_mean_vars=[] #20 amino acid composition+mean/var of 49 properties
    for i in range(len(lines)):
        if len(lines[i])>0 and lines[i][0]=='>':
            proseq=lines[i+1]
            #compute the amino acid composition
            aacontent=cvAA(proseq)
            #compute the k-property average,var
            proper49=mean_var_49(proseq)
            AACV_mean_vars.append(aacontent+proper49)
    # one id one list (20+49*2=118) 
    # 2d list
    return(AACV_mean_vars)




#if __name__=="__main__":
#    fr_n=sys.argv[1]
#    aa=api_aacv_pp49(fr_n)
#    print(aa)


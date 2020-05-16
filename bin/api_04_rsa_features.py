#!/usr/bin/env python
"""
this script is used to compute the featues of RSA
BE:Bd/Ed buired /Exposed
content of Bd, Ed
average of RSA in windows size 5,7,9,11,13,15,17,19,21
"""
import math,sys
from my_average import my_average

    
def get_rsa_seq(rsaValues,CUTOFF=0.25):
  "return rsa sequence"
  rsaseq=""
  for i in range(len(rsaValues)):
    if rsaValues[i]>CUTOFF:
      rsaseq=rsaseq+"1"
    else:
      rsaseq=rsaseq+"0"
  return rsaseq
#------------
def content_entropy(oneseq,flag="ss"):
  """input oneseq, return content """
  if flag=="ss":
    db={'C':0,'H':1,'E':2}
    f_count=[0 for i in range(3)]
    for i in range(len(oneseq)):
      f_idx=db[oneseq[i]]
      f_count[f_idx]=f_count[f_idx]+1
  elif flag=="rsa":
    db={'0':0,'1':1}
    f_count=[0,0]
    for i in range(len(oneseq)):
      if oneseq[i]=="0":
        f_count[0]=f_count[0]+1  # buried
      else:
        f_count[1]=f_count[1]+1  #exposed
  else:
    print( "no such flag")
    return 0.0
  #compute the content
  f_result=[]

  for i in range(len(db)): #compute the fraction
    f_p=f_count[i]*1.0/len(oneseq)
    f_result.append(f_p)
  return f_result
#-CHE+entropy
#Bd/Ed+entropy
#-----------------------
def combine_feature(seq1,seq2,flag1,flag2):
  """comput the compostion of AA in buried/exposed"""
  if flag1=="rsa"   and flag2=="seq" :
  # average of RSA in C, H, E
    #db={'C':0,'H':1,'E':2}
    AAdb={'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,\
    'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19,'U':20,'X':20}
    Bd_count=[0 for i in range(20)]
    Ed_count=[0 for i in range(20)]
    Bd_cv=[]
    Ed_cv=[]
    for i in range(len(seq2)):
      idx=AAdb[seq2[i]]
      if idx<20:
        if seq[i]=="0": #buried
          Bd_count[idx]=Bd_count[idx]+1
        elif seq[i]=="1": # exposed
          Ed_count[idx]=Ed_count[idx]+1
      else:
        print("### find X/U")
    Bdlen=sum(Bd_count)
    Edlen=sum(Ed_count)
    #compute the buried composition
    if Bdlen >0:
      for item in Bd_count:
        Bd_cv.append(item*1.0/Bdlen)
    else:# non-Buried residues
      for k in range(20):
        Bd_cv.append(0)
    #compute the exposed composition
    if Edlen >0:
      for item in Ed_count:
        Ed_cv.append(item*1.0/Edlen)
    else:# non-Buried residues
      for k in range(20):
        Ed_cv.append(0)
  return (Bd_cv,Ed_cv)

def sliding_windows(rsaValues,win_width):
  """ compute the average of sliding windows; return min,max"""  
  rsaseq=rsaValues
  win_num=len(rsaseq)-win_width+1
  if win_num >0:
    #set the intial number
    max_score=my_average(rsaseq[0:win_width])
    min_score=max_score
    for i in range( win_num ):
      onewin=rsaseq[i:i+win_width]
      oneave=my_average(onewin)
      # find the maximum
      if oneave > max_score:
        max_score=oneave
      if oneave < min_score:
        min_score=oneave
    #####
      #print(min_score,max_score,onewin)
  else: # win_num =0; or win_width larger than rsaValues
    max_score=0
    min_score=0
    print("### rsaValues len: %d < win-width %d" %(len(rsaValues),len(win_width)))
  return  (min_score,max_score)
# max_rsa_value min_rsa_values  average_raap_max_rsa, average_raap_min_rsa
#----------------2010-02-06--------

def api_rsa_features(seq,sss,xrsaValues):
  "input the rsa.dat ouput the rsa.csv"
  cutoffs=[0.25,0.75]
  win_width=[4,6,8,10,12,14,16,18,20,22]
  rsa_features=[]
  #rsa values ,change into 5digit to same as the query sequence 
  rsaValues=[]
  for item in xrsaValues:
    rsaValues.append(round(item,5))
  # different cutoff
  for onecutoff in cutoffs:
    rsa=get_rsa_seq(rsaValues,onecutoff)
    rsa_cv=content_entropy(rsa,"rsa")
    for onefeature in rsa_cv:
      # Bd/Ed_0.25  Bd/Ed_0.75 (4 values)
      rsa_features.append(onefeature)
      
  for j in range(len(win_width)):
    tmp_min,tmp_max=sliding_windows(rsaValues,win_width[j])
    # min/max_4/6/8/10/12/14/16/18/20/22; (20)
    rsa_features.append(tmp_min)
    rsa_features.append(tmp_max)

  return(rsa_features)

if __name__=="__main__":
  if len(sys.argv) <3:
    print("Usage rsa.features.py  fr_n  fw_n")
  else:
    fr_n=sys.argv[1]
    fw_n=sys.argv[2]
    fr=open(fr_n,'r')
    lines=[line.strip() for line in fr.readlines()]
    fr.close()
    #bline=[]
    seq=""
    sss=""
    for k in range(len(lines)):
      aline=lines[k].split(",")
      bline=[]
      for item in aline:
        if len(item)>0:
          bline.append(round(float(item),5))
      aa=api_rsa_features(seq,sss,bline) #   
      #print(aa) 

#"""

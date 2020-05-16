#!/usr/bin/env python
"""
this script is used to compute the featues of COBEprod
SS:C/H/E Coil,Helix,Strand
max segment
content of C, H, E.
# number of segements C,H,E.
Entropy of C/H/E

BE:Bd/Ed buired /Exposed
content of Bd, Ed
entropy of Bd, Ed.
average of RSA in windows size 5,6....17,18
Raap:
average of Raap  in widows size 5,6,...,17,18

combine the featuer SS,BE
average RSA in CHE
average Raap in CHE/BdEd
average Raap in B_CHE, Ed_CHE
---------2010-03-18-----------v.2
number of C_Bd H_Bd E_Bd
number of C_Ed H_Ed E_Ed
"""
import math,sys
#from Similarity import similarity_measure


def content_entropy(oneseq,flag="ss"):
  if flag=="ss":
    db={'C':0,'H':1,'E':2}
    f_count=[0 for i in range(3)]
    for i in range(len(oneseq)):
      f_idx=db[oneseq[i]]
      f_count[f_idx]=f_count[f_idx]+1
  elif flag=="rsa":
    db={'B':0,'E':1}
    f_count=[0,0]
    for i in range(len(oneseq)):
      if float(oneseq[i])<0.25:
        f_count[0]=f_count[0]+1  # buried
      else:
        f_count[1]=f_count[1]+1  #exposed
  else:
    print( "no such flag")
    return 0.0
  entropy=0
  f_result=[]

  for i in range(len(db)): #compute the fraction
    f_p=f_count[i]*1.0/len(oneseq)
    f_result.append(f_p)
    if(f_p>0):
      entropy=entropy+(-1)*f_p*math.log(f_p,2)
    else:
      entropy=entropy+0
  f_result.append(entropy)

  return f_result
#-CHE+entropy
#Bd/Ed+entropy

def average_rsa_be(rsaseq):
  f_result=[]
  bvalues=[]
  evalues=[]
  for i in range(len(rsaseq)):
    if rsaseq[i]<0.25:
      bvalues.append(rsaseq[i])
    else:
      evalues.append(rsaseq[i])
  [blen,elen]=[len(bvalues),len(evalues)]
  if blen>0:
    f_result.append(sum(bvalues)*1.0/blen)
  else:
    f_result.append(0.0)
  if elen>0:
    f_result.append(sum(evalues)*1.0/elen)
  else:
    f_result.append(0.0)
  return f_result
#---averag RSA of Buried , Exposed.
#rsaseq_tmp="0.760,0.591,0.414,0.448,0.326,0.511,0.210,0.505,0.511,0.301,0.313,0.384,0.260,0.294,0.431,0.384,0.262,0.428,0.616,0.657,".split(",")
#rsaseq=[float(rsaseq_tmp[one]) for one in range(len(rsaseq_tmp)-1)]
#print average_rsa(rsaseq)
#proseq="EHHESTWSDAYPYSKRMAEK"
#------------
def ss_segment(oneseq):
  segments=[]
  mychar=oneseq[0]
  myseg=""
  for i in range(len(oneseq)):
    if oneseq[i]==mychar:
      myseg=myseg+mychar
      #myse.append(oneseq[i])
    else:
      mychar=oneseq[i]
      segments.append(myseg)
      myseg=""
      myseg=myseg+mychar
      #myseg.append(mychar)
    if i==len(oneseq)-1:
      segments.append(myseg)
  return segments
#--segments of one sequence
#a1=ss_segment("ABBBCCEEEEEXXO")
#print a1
#----------------------
def ss_nums(mysegments):
  SSdb={'C':0,'H':1,'E':2}
  f_result=[0,0,0]
  for item in  mysegments:
    idx=SSdb[item[0]]
    f_result[idx]+=1
  return f_result
#------------------
#a1=ss_segment("CCCHHEEEEEC")
#print a1
#print ss_nums(a1)
#-----------------------
def combine_feature(seq1,seq2,flag1,flag2):
  if flag1=="rsa"   and flag2=="ss" :
  # average of RSA in C, H, E
    db={'C':0,'H':1,'E':2}
    mycount=[[] for i in range(3)]
    myresult=[]
    for i in range(len(seq2)):
      mycount[db[seq2[i]]].append(float(seq1[i]))
    for i in range(3):
      if len(mycount[i])>0:
        myresult.append(sum(mycount[i])/len(mycount[i]))
      else:
        myresult.append(0.0)

  return myresult
# rsa_C, rsa_H, rsa_E
#seq2="CCECCCCCCCCCCCCCCCCC"
#seq1="0.536,0.478,0.215,0.463,0.466,0.586,0.495,0.623,0.571,0.466,0.680,0.496,0.358,0.546,0.343,0.466,0.557,0.520,0.563,0.742,".split(",")
#print combine_feature(seq1,seq2,"rsa","ss")
#-------------------
def average_Raap(Raap,seq):
  AAdb={'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,\
    'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19}
  #Raap=scale
  if len(seq)<2:
    return 0.0
  else:
    total=0.0
    #f_values=[]
    for i in range (len(seq)-1):
      idx=AAdb[seq[i]]*20+AAdb[seq[i+1]]
      total+=Raap[idx]
    myave=total*1.0/(len(seq)-1)
  return myave

#map the sequence to the average of Raap values
#proseq="VRWRRKSSDRKGGSYSQAAS"
#average_Raap
#------------
def sd_Raap(Raap,seq,mean_value):
  AAdb={'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,\
    'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19}
  #Raap=scale
  if len(seq)<2:
    return 0.0
  else:
    total=0.0
    for i in range (len(seq)-1):
      idx=AAdb[seq[i]]*20+AAdb[seq[i+1]]
      total+=(Raap[idx]-mean_value)*(Raap[idx]-mean_value)
    return math.sqrt(total*1.0/(len(seq)-1))
#-------------------------- 
def sliding_window_2D(Scale,seq,win_width):
  AA={'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,\
    'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19}
  #Raap=scale
  if len(seq)< 2:
    return 0.0
  f_result=[]
  max_score=-20000
  min_score=10000
  win_num=len(seq)-win_width+1
  for i in range( win_num ): # win_num: number of windows
    onewin=seq[i:i+win_width]
    oneave=average_Raap(Scale,onewin)
    #print oneave
    if oneave > max_score:
      max_score=oneave
      #max_seq=onewin
    if oneave < min_score:
      min_score=oneave
      #min_seq=onewin
  f_result.append(max_score)
  
  f_result.append(min_score)
  
  return  f_result       # average of the Raap value in all windows
#----f_result[0]: max_Raap_values min_Raap_values
#
def sliding_window_1D(Raap,proseq,rsaseq,win_width):
  AA={'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,\
    'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19}
  #seq=rsaseq
  f_result=[]
  max_score=-10000
  min_score=10000
  win_num=len(rsaseq)-win_width+1
  for i in range( win_num ):
    onewin=rsaseq[i:i+win_width]
    #print onewin
    oneave=sum(onewin)*1.0/win_width
    #print oneave
    if oneave > max_score:
      max_score=oneave
      max_star=i
    if oneave < min_score:
      min_score=oneave
      min_star=i
  max_seq=proseq[max_star:max_star+win_width]
  min_seq=proseq[min_star:min_star+win_width]
  average_raap_max_rsa=average_Raap(Raap,max_seq)
  average_raap_min_rsa=average_Raap(Raap,min_seq)
  f_result.append(max_score)
  f_result.append(min_score)
  f_result.append(average_raap_max_rsa)
  f_result.append(average_raap_min_rsa)
  return  f_result
# max_rsa_value min_rsa_values  average_raap_max_rsa, average_raap_min_rsa
#----------------2010-02-06--------
#------max-segments-----
def max_segment( Raap, mysegments,proseq,ssseq,rsaseq,conseq): #2010-03-30
  AAdb={'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,\
    'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19}

  my_results=[]
  SSdb={'C':0,'H':1,'E':2}
  ssmax=[[0,""],[0,""],[0,""]]
  for i in range( len(mysegments)):
    item=mysegments[i]
    myidx=SSdb[item[0]]
    if len(item) > ssmax[myidx][0]:
      ssmax[myidx][0]=len(item) #record the length
      ssmax[myidx][1]=item     # record the sequence
  #print ssmax
  # find the correspond the sequence
  for i in range(3):
    len_max=ssmax[i][0] # segment i length. i=C,H,E.\
    #print len_max
    rsa_value=0.0
    con_value=0.0 #2010-03-30
    raap_value=0.0
    if len_max > 1:
      start=ssseq.find(ssmax[i][1])
      if start >0 : # find it
        proseg=proseq[start:start+len_max]
        rsaseg=rsaseq[start:start+len_max]
        conseg=conseq[start:start+len_max] # 2010-03-30
        rsa_value=sum(rsaseg)*1.0/len_max
        con_value=sum(conseg)*1.0/len_max # 2010-03-30
        Raap_items=[]
        for k in range(len_max-2+1):
          myidx=AAdb[proseg[k]]*20+AAdb[proseg[k+1]]
          Raap_items.append(Raap[myidx])
        raap_value=sum(Raap_items)/(len_max-1)
    elif len_max==1:
      start=ssseq.find(ssmax[i][1])
      rsa_value=rsaseq[start]
      con_value=conseq[start] # 2010-03-30
      raap_value=0.0
    my_results.append(rsa_value)# add the average of rsa
    my_results.append(con_value) # 2010-03-30
    my_results.append(raap_value)# add the average of Raap

  return my_results
# C_rsa,C_raap,H_rsa,H_raap,E_rsa,E_raap
#----------------------------------
def Raap_feature(Raap, proseq,ssseq,rsaseq):
  myresults=[]
  SSdb={'C':0,'H':1,'E':2}
  AAdb={'A':0,'C':1,'D':2,'E':3,'F':4,'G':5,'H':6,'I':7,'K':8,'L':9,\
    'M':10,'N':11,'P':12,'Q':13,'R':14,'S':15,'T':16,'V':17,'W':18,'Y':19}
  SSpro_seg=["","",""]
  BEpro_seg=["",""]
  SSBpro_seg=["","",""] #buried
  SSEpro_seg=["","",""] #exposed

  for i in range(len(ssseq)):
    idx=SSdb[ssseq[i]]
    #SSpro_seg[idx].append(proseq[i])
    SSpro_seg[idx]+=proseq[i]
    if rsaseq[i]<0.25: #buried
      #BEpro_seg[0].append(proseq[i])
      #SSBpro_seg[idx].append(proseq[i])
      BEpro_seg[0]+=proseq[i]
      SSBpro_seg[idx]+=proseq[i]
    else:
      #BEpro_seg[1].append(proseq[i])
      #SSEpro_seg[idx].append(proseq[i])
      BEpro_seg[1]+=proseq[i]
      SSEpro_seg[idx]+=proseq[i]
  #print SSpro_seg
  for i in range(3):
    #myresults.append(sliding_window_2D(Raap,SSpro_seg[i],len(SSpro_seg[i]) )) # average Raap in C,H,E
    myresults.append(average_Raap(Raap,SSpro_seg[i]))
  for i in range(2): # average Raap in Buried Exposed
    #myresults.append(sliding_window_2D(Raap,BEpro_seg[i],len(BEpro_seg[i]) ))
    myresults.append(average_Raap(Raap,BEpro_seg[i]))
  for i in range(3): #average Raap in C,H,E in Buried
    #myresults.append(sliding_window_2D(Raap,SSBpro_seg[i],len(SSBpro_seg[i]) ))
    myresults.append(average_Raap(Raap,SSBpro_seg[i]))
  for i in range(3): # average Raap in C,H,E in Exposed
    #myresults.append(sliding_window_2D(Raap,SSEpro_seg[i],len(SSEpro_seg[i]) ))
    myresults.append(average_Raap(Raap,SSEpro_seg[i]))

  return myresults
#--------------------------------
def num_features(ssseq,rsaseq):
  SSdb={'C':0,'H':1,'E':2}
  myresult=[]
  ssnum=[0,0,0]
  benum=[0,0]
  bessnum=[[0,0,0],[0,0,0]]
  for item in ssseq:
    ssnum[SSdb[item]]+=1
  for i in range(len(rsaseq)):
    if rsaseq[i] < 0.25:
      benum[0]+=1
      idx=0
      bessnum[0][SSdb[ssseq[i]]]+=1
    else:
      benum[1]+=1
      idx=1
      bessnum[1][SSdb[ssseq[i]]]+=1
  #write 
  #for i in range( len(ssnum)):
  #  myresult.append(ssnum[i])
  #for i in range( len(benum)):
  #  myresult.append(benum[i])
  for i in range(2):
    for j in range(3):
      myresult.append(bessnum[i][j])
  return myresult
# number of C, H, E
# number of Bd Ed ; C/H/E_Bd; C/H/E_Ed

#-----------Main---Function----------
input_f=sys.argv[1] #infut file ../structure_results/epitopia.stru.test.ss.rsa.rh.txt
output_f=sys.argv[2] #output file #../structure_results/epitopia.stru.test.no.similar.csv
label=sys.argv[3]
#input_f="../data/non.ss.rsa.txt"#"../data/vol.ss.rsa.txt"#"vol.ss.rsa.txt" #remembe to change label
#output_f="../data/non.ss.rsa.csv"#"../data/vol.features.csv"##"vol.features.csv"
#Note: rember change label in Line367
f1=open(input_f,'r')
#f1=open("bcpred.n.ss.rsa.txt",'r') ../data/bcpred20.train.ss.rsa.rh.txt
data=[line.strip() for line in f1.readlines()]
f1.close()
fw=open(output_f,'w')
##fw=open("bcpred.n.svm",'w') ../data/bcpred20.train.no.similar.csv
#label="1"

#u"""
datasize=len(data)

SSdic={'C':0,'H':1,'E':2}

for i in range(datasize):
  if len(data[i])>0 and data[i][0]=='>':
    #label="1" 
    #if data[i].split(".")[1]=="p":
    #  label="+1"
    #else:
    #  label="-1"
    seq=data[i+1]
    sss=data[i+2]
    trsa=data[i+3].split(",") # remove last item
    #tcon=data[i+4].split(",") #conservation score
    #rsa sequence
    rsa=[float(trsa[item]) for item in range (len(trsa)-1)]
    #conx=[float(tcon[item]) for item in range (len(tcon)-1)]
    # ss fetures: fraction of CHE, entropy (3+1)
    f1_ss_cv_entropy=content_entropy(sss,"ss")
    # Be features: fraction of Bd, Ed, entropy (2+1)
    f2_rsa_cv_entropy=content_entropy(rsa,"rsa")
    # 
    #avaerg of RSA of buried and Exposed (2)
    f3_ave_rsa_be=average_rsa_be(rsa)
    #f31_ave_con_be=sum(conx)/len(conx)  # 2010-03-30 conx means conversation
    #print "f31 ok"
    #ss segment number  (1)
    mysegments=ss_segment(sss)
    f4_ss_seg_num=len(mysegments)
    # number of segments in C/H/E  (3)
    f5_ss_num_CHE=ss_nums(mysegments)
    #average of RSA of C/H/E (3)
    f6_rsa_ss=combine_feature(rsa,sss,"rsa","ss")
    
    # number of C/H/E_Bd  C/H/E_Ed (6)
    f13_nums=num_features(sss,rsa)    
    #fw.write("%s," %label)  ########################## 2012-04-23 just for train
    #idx_svm=1
	#SS: 1,2,3,4
    for k in range(len(f1_ss_cv_entropy)):
      fw.write("%.5f," %(f1_ss_cv_entropy[k]))
	#RSA 5,6,7
    for k in range(len(f2_rsa_cv_entropy)):
      fw.write("%.5f," %(f2_rsa_cv_entropy[k]))
	#RSA: 8,9
    for k in range(len(f3_ave_rsa_be)):
      fw.write("%.5f," %f3_ave_rsa_be[k])
    
    #2010-03-30
    #fw.write("%.3f," %f31_ave_con_be)
    #SS:10
    fw.write("%d," %(f4_ss_seg_num))
	#SS:11,12,13
    for k in range(len(f5_ss_num_CHE)):
      fw.write("%.5f," %f5_ss_num_CHE[k])
	#RSA:14,15,16
    for k in range( len(f6_rsa_ss)):
      fw.write("%.5f," %f6_rsa_ss[k])
	#SS:17,18,19; 20,21,22
    for k in range(len(f13_nums)):
      fw.write("%d," %(f13_nums[k]))
    fw.write("\n")


#"""

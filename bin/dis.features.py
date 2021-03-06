#!/usr/bin/env python
u"""
this script is used to compute the featues based on 
SS:C/H/E Coil,Helix,Strand
"""
import math,sys
#from Similarity import similarity_measure
#------------
def ss_segment(oneseq):
  #return the ss-segments		
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
def total_ss(segments):
  return len(segments)
#print a1
#----------------------
def ss_nums(mysegments):
  # number of ./* fragments .means order * disorder
  SSdb={'.':0,'*':1}
  f_result=[0,0]
  for item in  mysegments:
    idx=SSdb[item[0]]
    f_result[idx]+=1
  return f_result
#-----------------
def percent_min_max_ss(mysegments,proseq):
  # percent of min /max  C/H/E fragments. i.e. max/min fragments divided by sequence length
  SSdb={'.':0,'*':1}
  seqlen=len(proseq)
  f_result=[[0,0],[0,0]] #min/max for C/H/E
  for item in mysegments:
    seglen=len(item) # fragment/segments length
    idx=SSdb[item[0]]
    if f_result[idx][0] > seglen: #min < seglen
      f_result[idx][0]=seglen
    if f_result[idx][1] < seglen: #max < seglen
      f_result[idx][0]=seglen
  percent_f=[]
  for i in range(len(f_result)):
    percent_f.append(f_result[i][0]*1.0/seqlen)
    percent_f.append(f_result[i][1]*1.0/seqlen)
   
  return percent_f

#----------------
def cv_ss(oneseq):
  #return the  percentage of C/H/E
  SSdb={'.':0,'*':1}
  cvss=[0,0]
  results=[]
  seqlen=len(oneseq)
  if seqlen!=0:
    for i in range(seqlen):
      idx=SSdb[oneseq[i].upper()]
      cvss[idx]=cvss[idx]+1
    for i in range(len(cvss)):
      results.append(cvss[i]*1.0/seqlen)
  else:
    print("empty sequence")
    results=[0,0]
  return results
#-------------------

# max_rsa_value min_rsa_values  average_raap_max_rsa, average_raap_min_rsa
#----------------2010-02-06--------
#-----------Main---Function----------
def ss_bash(input_f, output_f):

  f1=open(input_f,'r') # ss/dis-order fasta 
  data=[line.strip() for line in f1.readlines()]
  f1.close()
  fw=open(output_f,'w')
  
  datasize=len(data)
  SSdic={'.':0,'*':1}
  ss_features=[]
  for i in range(datasize):
    if len(data[i])>0 and data[i][0]=='>':
      seq=data[i+1]
      sss=data[i+2]
      #ss segements CC HHH EEEE CC
      mysegments=ss_segment(sss)

      # total number of ss fragments (1)
      ss_seg_num=total_ss(mysegments)
      fw.write("%d," %ss_seg_num)

      # number of segments in C/H/E  (2)
      ss_num_CHE=ss_nums(mysegments)
      for k in range(len(ss_num_CHE)):
        fw.write("%d," %(ss_num_CHE[k]))
      #min/max segment fraction of order/ Disordered  (2*2=4 )
      minmax_ss=percent_min_max_ss(mysegments,seq)
      for k in range(len(minmax_ss)):
        fw.write("%.5f," %minmax_ss[k])
      
      # composition of ss (2)
      percent_ss=cv_ss(sss)
      for k in range(len(percent_ss)):
        fw.write("%.5f," %percent_ss[k])
      fw.write("\n")

if __name__=="__main__":
  if len(sys.argv)<3:
    print("Usage dis.features.py fr_n fw_n\nfr_n:id;seq;dis;")
  else:
    fr_n=sys.argv[1]
    fr_w=sys.argv[2]
    ss_bash(fr_n,fr_w)

#"""

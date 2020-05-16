# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 17:45:23 2019
### add constrains 
@author: Cliff-XPS
"""
import os
#import numpy as np
import sys
id=sys.argv[1]
type=sys.argv[2]
outputdir='./output/%s'%id

aname=['ionchannel-ligand-potassium','ionchannel-ligand-calcium',
           'ionchannel-ligand-sodium','ionchannel-ligand-anion',
           'ionchannel-voltage-potassium','ionchannel-voltage-calcium',
           'ionchannel-voltage-sodium','ionchannel-voltage-anion']   
pred_out=['ligand-gated channel for potassium','ligand-gated channel for calcium',
           'ligand-gated channel for sodium','ligand-gated channel for anion',
           'voltage-gated channel for potassium','voltage-gated channel for calcium',
           'voltage-gated channel for sodium','voltage-gated channel for anion']    
dict_out=dict(zip(aname,pred_out))

def read_3predictions():
    ### read in probabilities from PSIONplus three models output file
    models=["ION","VLG","VGS"]
    probs=[] #non-ion;ion;lig;vol;K,Ca;Na;Anion #multi-label order
    #psionplus order:k, anion,ca,na
    for amodel in models:
        afrn=outputdir+"/output.%s.psionplus.predict" %amodel  #change it in PSIONplus4*.sh the same time 
        if os.path.isfile(afrn):
            fr=open(afrn,'r')
            for line in fr.readlines():
                if line[0]!='#':
                    if amodel !='VGS':
                        apred,negP,posP=line.split()
                        probs.append(float(negP))
                        probs.append(float(posP))
                    else:
                        apred,kp,anionp,cap,nap=line.split()
                        probs.append(float(kp))
                        probs.append(float(cap))
                        probs.append(float(nap))
                        probs.append(float(anionp)) #order as multi-label order
        else:
            print("Can not find %s" %afrn)
    #non-ion;ion;lig;vol;K,Ca;Na;Anion #multi-label order        
    return(probs)
#probs=read_3predictions()    
#print(probs)

def my_joint_probs(probs):
    #### compute the ionchannel-lig/vol-k/ca/na/aion probabilities
    #### non-ion,ion,ligand,voltage,4types
    jointProbs=[]
    nonP=probs[0] #probability of non-ion channel
    ionP=probs[1] #probability of ion channel
    ligP=probs[2]
    volP=probs[3]
    typP=probs[4:]
    if ionP > nonP :
        #only selected ionchannel-lig/vol-k/ca/na/anion
        for k in range(len(typP)):
            jointProbs.append(ligP*typP[k])
        for k in range(len(typP)):
            jointProbs.append(volP*typP[k])
    else:
        #non-ion-lig/vol-k/ca/na/anion
        pass
    ### lig/vol-k/ca/na/anion
    return(jointProbs)
def my_constrains(alist,aquery):
    addFlag=True
    #if query is  
    ## do not add the extra label 
    ## if both atom type and channel type are different.
    #AtomTypes=[] # K/Ca/Na/Cl
    #ChannelTypes=[] #ligand/types
    if len(alist)>0:
        qatom=aquery.split('-')[1].strip()
        qchan=aquery.split('-')[2].strip()
        
        for item in alist:
            atom_type=item.split('-')[1].strip()
            channel_type=item.split('-')[2].strip()
            if (qatom!=atom_type) and (qchan!=channel_type):
                addFlag=False
                break
    else:
        addFlag=True
    return(addFlag)
def gen_pred_above_cutoff(xprobs,cutoff):
    ## nonIon,Ion,Lig,Vol,1':'potassium','2':'calcium','3':'sodium','4':'anion'
    jointProbs=my_joint_probs(xprobs)
    ### jointProbs: lig-k/ca/na/cl; vol-k/ca/na/cl/
    preds=[]
    for idx,aprob in enumerate(jointProbs):
        if aprob>cutoff:
            apred=aname[idx]
            ### 2019-10-23 add constrains
            addFlag=my_constrains(preds,apred) 
            if addFlag:
                preds.append(apred)
    ## if there is no score >=0 cutoff, choose max probability
    if len(preds)==0:
        maxidx=np.argsort(jointProbs)[::-1][0] #decrease
        preds.append(aname[maxidx])
    return(preds)
    
def gen_pred_single(xprobs):
    ### return xprobs >= cutoff
    ## nonIon,Ion,Lig,Vol,1':'potassium','2':'calcium','3':'sodium','4':'anion'
    jointProbs=my_joint_probs(xprobs)
    ### jointProbs: lig-k/ca/na/cl; vol-k/ca/na/cl/
    
    #maxIdx=np.argsort(jointProbs)[::-1] #revser []
    #topIdx=maxIdx[0:topN]
    idx=jointProbs.index(max(jointProbs))
    apred=aname[idx]
    ## if there is no score >=0 cutoff, choose max probability
    return(dict_out[apred])
    

def  get_sequential_pred(probs,cutoff):
    # get the sequenctial prediction
    seq_pred="" #mutlti-label; score >cutoff
    nonP=probs[0] #probability of non-ion channel
    ionP=probs[1] #probability of ion channel
    ligP=probs[2]
    volP=probs[3]
    typP=probs[4:]#K,Ca;Na;Anion #multi-label order

    if nonP>=ionP:
        seq_pred="non_ion"
    else:
        multi_labels0=gen_pred_above_cutoff(probs,cutoff) #20191108
        multi_labels=[dict_out[k] for k in multi_labels0]
        seq_pred=";".join(multi_labels)

    return(seq_pred)
  
          
def main():
    ### get the squential prediction above cutoff
    probs=read_3predictions()
    if type=='5':
        cutoff=0  #40%  top-40% scores
        seq_pred=get_sequential_pred(probs,cutoff)
        fwn=outputdir+'/output.SEQ.psionplusM.single_label.predict'              
        fw=open(fwn,'w')
        fw.write("%s\n" %seq_pred)
        fw.close()
    elif type=='4':
        pred=gen_pred_single(probs)
        fwn=outputdir+'/output.SEQ.psionplusM.multi_label.predict'              
        fw=open(fwn,'w')
        fw.write("%s\n" %pred)
        fw.close()
    print("### output is %s" %fwn)
    return(0)
if __name__ == "__main__":
	main()



  

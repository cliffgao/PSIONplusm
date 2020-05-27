#!/usr/bin/env python
import os, sys, gc, getopt
import multiprocessing
import time
import argparse
from argparse import RawTextHelpFormatter

############################
pwd=os.getcwd()
#DISOPREDdir=pwd+'/software/DISOPRED'
#SPINEXdir=pwd+'/software/spineXpublic'
#PSIPREDdir=pwd+'/software/psipred'
DISOPREDdir=pwd+'/lib/DISOPRED'
SPINEXdir=pwd+'/lib/spineXpublic'
PSIPREDdir=pwd+'/lib/psipred'
bindir=pwd+'/bin'
NCBIdir="/home/yangjy/I-TASSER/I-TASSER4.4/blast/bin"
db="/library/nr/nr"
##############################

    
    


def diso(outputdir):
    if not os.path.exists('%s/target.diso'%(outputdir)):
        os.system('%s/run_disopred.pl %s/target.fasta %s %s'%(DISOPREDdir,outputdir,NCBIdir,db))

   
def spX(outputdir):
    if not os.path.exists('%s/spXout/target.spXout'%(outputdir)):
        f=open('%s/list'%outputdir,'w')
        f.write('target')
        f.close()
        os.system('cp target.fasta target')
        os.system('%s/spX.pl list %s %s %s'%(SPINEXdir,outputdir,NCBIdir,db))
def ss(outputdir):
    if not os.path.exists('%s/target.ss'%(outputdir)):
        os.system('%s/runpsipred target.fasta %s %s'%(PSIPREDdir,db,NCBIdir))
###
def get_features(outputdir):
    os.chdir('%s'%outputdir)
    p1 = multiprocessing.Process(target=diso(outputdir))
    p2 = multiprocessing.Process(target=spX(outputdir))
    p3 = multiprocessing.Process(target=ss(outputdir))
    p1.start()
    p2.start()
    p3.start()
    p1.join()
    p2.join()
    p3.join()



if __name__=="__main__":
    if len(sys.argv)<2:
        print ('usage: python pred.py [-h] [-i fastafile] [-m model]')
        sys.exit(0)
    helptxt='''
the model used to predict.
options:
1 PSIONplus:  discrimination for ion channel and non-ion channel (ION)
2 PSIONplus:  discrimination for voltage-gated and ligand-gated channel (VLG)
3 PSIONplus:  discrimination for four types ion channel (VGS)
4 PSIONplusm: sequential prediction for single-label prediction
5 PSIONplusm: sequential prediciton for multi-label prediction)
'''
    parser = argparse.ArgumentParser(description='Predict ion channels and their types from protein sequence.',
            formatter_class=RawTextHelpFormatter)    
    parser.add_argument("-i","--infile", help="the protein sequence in fasta format.")
    parser.add_argument("-m","--model",  choices=['1','2','3','4','5'], help=helptxt)
    args = parser.parse_args()        
    aid=args.infile
    atype=args.model
    ##
    if '.' in aid:
        aid=aid.split('.')[0]
    else:
        pass
    outputdir=pwd+'/output/%s'% aid
    if not os.path.exists(outputdir):
        os.mkdir(outputdir)
    else:
        #print("%s exists" %outputdir)
        pass
    if not os.path.isfile(args.infile):
        print("Cant not find %s" %args.infile)
        sys.exit(0)
    # copy the fasta sequence in the outputdir as target.fasta 
    os.system("cp %s  %s/target.fasta" %(args.infile,outputdir))
    get_features((outputdir))  
    if os.path.exists('%s/target.mat'%(outputdir)) and not os.path.exists('%s/target.pssm'%(outputdir)):
        os.system('mv target.mat target.pssm')
    os.chdir(pwd)
    ## 
    #os.system("chmod +x ./bin/PSIONplus4multi_labels.sh")
    os.system("chmod +x ./bin/*")
    if atype in ['1','2','3']:  
        os.system('%s/bin/PSIONplus4multi_labels.sh %s %s'%(pwd,aid,atype))
        # output the message 
        mdlsdb={'1':'ION','2':'VLG','3':'VGS'}
        try:
            xmdl=mdlsdb[atype]
        except KeyError:
            print("error model type(1/2/3)")
            xmdl="X"
        currentdir="./output/%s" %aid
        print("### output is %s/output.%s.psionplus.predict" %(currentdir,xmdl))
    elif atype in ['4','5']: 
        os.system('%s/bin/PSIONplus4multi_labels.sh %s 1'%(pwd,aid))
        os.system('%s/bin/PSIONplus4multi_labels.sh %s 2'%(pwd,aid))
        os.system('%s/bin/PSIONplus4multi_labels.sh %s 3'%(pwd,aid))
        os.system('python %s/bin/wb_get_sequential_pred_with_probs.py %s %s'%(pwd,aid,atype))



    

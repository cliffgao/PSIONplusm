#!/usr/bin/env python

"""scale the svm model in to  -1 1 for test dataset """
import sys
def get_scale_file(fr_n):
    #input the fr_n, 1st min; 2nd max;
    fr=open(fr_n,'r')
    lines=[line.strip() for line in fr.readlines()]
    fr.close()
    mymin=[]
    mymax=[]
    #----------- new version
    for i in range(len(lines)):
        if lines[i]:
            iteMin,iteMax=lines[i].split()
            mymin.append(float(iteMin))
            mymax.append(float(iteMax))
    return(mymin,mymax)

def scale_model(filename,fw_n,item_mi,item_ma):
    fr=open(filename,'r')
    lines=[line.strip() for line in fr.readlines()]
    fr.close()
    fw=open(fw_n,'w')

    size=len(lines)
    features=[[] for i in range(size)]
    label=[]
    for i in range(size):
        if len(lines[i])>0 :
            eachline=lines[i].split()
            label.append(eachline[0])
            
            for j in range(1,len(eachline)):
                item=float(eachline[j].split(":")[1])
                features[i].append(item)
#    #scale the features
    for i in range(len(features)):
        fw.write("%s " %label[i])
        #print("%s " %label[i]),
        svm_idx=1
        for j in range(len(features[i])):
            if item_ma[j]==item_mi[j] :  # if the fenmu is zero
                item=item_ma[j]
            else:
                item=(features[i][j]*2.0-(item_mi[j]+item_ma[j]))/(item_ma[j]-item_mi[j])
            #print("%d:%.5f " %(svm_idx,item)),
            fw.write("%d:%.4f " %(svm_idx,item))

            svm_idx+=1
        fw.write("\n")
        #print("")
    # features, min,max
    return (item_mi,item_ma)

if __name__=="__main__":
    if len(sys.argv)!=4:
        print("Usage my_scale_model.py  svm_fr_n min_max_fr_n(1st min/2nd max/sep by blank) fw_n")
    else:
        fr_name=sys.argv[1]
        scalefile=sys.argv[2]
        fw_n=sys.argv[3]
        mymin,mymax=get_scale_file(scalefile)
        scale_model(fr_name,fw_n,mymin,mymax)
        #print mymin
        #print mymax


                




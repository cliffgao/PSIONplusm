#!/bin/bash
#PSIONplus software method
id=$1
my_type=$2
datadir="./data/"
exedir="./bin/"
resultdir="./output/"${id}/
modeldir="./data/train.model/"
scaledir="./data/scale.id/"
#svmmodeldir="./data/libsvm/"
svmmodeldir="./lib/libsvm/"
blastdir="./data/blast.train/"
blastout="./data/blast.test/"
#blastalldir="./software/blast-2.2.26/bin/"
blastalldir="./lib/blast-2.2.26/bin/"

# input the fasta files and predictoin type
#1: for discrimination ion channel(1) and non-ion channel(-1)
#2: for discriminatoin voltage-gated ion channel(1) and ligand-gated ion channel(-1)
#3: for discrimination of four types of voltage-gated ion channel, i.e. potassium (K, 1), anion (Anion, 2), calcium( Ca, 3), sodium (Na, 4)
if [ $# -lt 2 ]
then
echo "Usage: PSIon.sh  fasta_file  types(1/2/3)"
echo "1: for discrimination for ion channel and non-ion channel (ION)"
echo "2: for discrimination for voltage-gated and ligand-gated channel (VLG)"
echo "3: for discrimination for voltage-gated ion channel subtypes (VGS)"
else

fasta_f=${resultdir}target.fasta
#"3" #"2" #"1"

if [ ${my_type} == "1" ]
then

	oneModel="ION" #"ion_non"

	e_value="0.001"
	p_cutoff="0.4"
	message="label: ion channel(1), non-ion channel (-1)"
elif [ ${my_type} == "2" ]
then
	oneModel="VLG" #"vol_lig"

	e_value="10"
	p_cutoff="0.5"
	message="label: voltage-gated ion channel(1), ligand-gated ion channel(-1)"
elif [ ${my_type} == "3" ]
then
	oneModel="VGS" #"types"

	e_value="0.001"
	p_cutoff="0.3"
	message="label:potassium (K, 1), anion (Anion, 2), calcium( Ca, 3), sodium (Na, 4)"
fi
### svm model files 
f_scale_info=${scaledir}train.${oneModel}.svm.info.scale
f_select_idx=${scaledir}train.${oneModel}.select.id

f_select_features_scale_model=${resultdir}test.${oneModel}.scale.select.svm  #query sequence model 
f_trained_model=${modeldir}train.${oneModel}.model  #PSION model 
train_fasta_file=${blastdir}Train.${oneModel}.txt    #training dataset for blast 
#ouput of blast hits
out_blast=${blastout}"out_blast"
out_blast_fasta=${blastout}"out_blast_fasta"
psion_out=${resultdir}zzz.output.${oneModel}.psion.predict  #prediction from just psion
#prob_txt=${resultdir}output.${oneModel}.combine.predict  #prediction from psionplus( combined blast)
prob_txt=${resultdir}output.${oneModel}.psionplus.predict  #prediction from psionplus( combined blast)
blast_out_flag=False #not ouput the blast prediction
#:<<BLOCK

${exedir}output_svm_mdl.py ${fasta_f} ${resultdir} ${f_scale_info}  ${f_select_idx}  ${f_select_features_scale_model}

# predict 
chmod +x ${svmmodeldir}svm-predict
${svmmodeldir}svm-predict -b 1 ${f_select_features_scale_model}  ${f_trained_model}  ${psion_out}




###--Running BLAST results

#running blast
#if [ -f "blastall" ]; then
chmod +x ${blastalldir}blastall
${blastalldir}blastall -p blastp -i ${fasta_f} -d ${train_fasta_file} -o ${out_blast}  -e ${e_value}
#else
#echo "Please install blastall/blastp"
#fi

#prepare blast into fasta format
${exedir}8get_hit_from_blast.py ${out_blast} ${out_blast_fasta} ${e_value}
#combine the probabilitiy

${exedir}9compute.prob.cutoff4multi_labels.py  ${out_blast_fasta} ${prob_txt} ${my_type} ${psion_out} ${e_value} ${blast_out_flag}
#remove the tmp file
if [ ${blast_out_flag} == "False" ]
then
rm ${resultdir}zzz.*
fi

#if [ ${my_type} == "1" -o ${my_type} == "2" -o ${my_type} == "3" ]
#then    
#    echo "### Output is" ${prob_txt}
#elif [ ${my_type} == "4" -o ${my_type} == "5" ]
#then 
#    echo "Sequential Prediction"
#fi
### if my_type is 4/5 donot output message 
#echo "### "${message}
fi

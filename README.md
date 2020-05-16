# PSIONplusm
This is a stand-alone version of PSIONplusm for accurate multi-label prediction of ion channels and their types.

# webserver
You can also use our server at:  
https://yanglab.nankai.edu.cn/PSIONplusm/


# The software need to install:
## 1.SPINE X  
Download the source code of Real SPINE X and install it.  
http://sparks-lab.org/SPINE-X/  
## 2.PSIPRED 3.3  
Download the PSIPRED from:  
http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/
## 3.DISOPRED  
Download the DISOPRED 3.3 from:  
http://bioinfadmin.cs.ucl.ac.uk/downloads/DISOPRED/  
## 4. PSIBLAST and BLAST  
Download the PSIBLAST from  
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/2.2.6/  

## Run PSIPRED,DISOPRED, SPINE X, PSIBLAST

1. run PSIPREd and get the result named target.diso.  Copy the target.diso into ./output/example/  
2. run PSIBLAST and get the result named target.pssm. Copy the target.pssm into ./output/example/  
3. run SPINE X and get the result named target.spXout. Copy the target.spXout into ./output/example/spXout  
4. run PSIPRED and get the result named target.ss. Copy the target.ss into ./output/example/  

# Change the Parameters

1.open the pred.py  
2.change the NCBIdir="/home/yangjy/I-TASSER/I-TASSER4.4/blast/bin" with your psiblast path (line 13)  
3.change the db="/library/nr/nr" with your nr database path (line14)  


# run the PSIONplusm  

Usage: pred.py  -i id  -t type(1/2/3/4/5)  

1 PSIONplus: discrimination for ion channel and non-ion channel (ION)  
2 PSIONplus: discrimination for voltage-gated and ligand-gated channel (VLG)  
3 PSIONplus: discrimination for four types ion channel (VGS)  
4 PSIONplusm: sequential prediction for single-label prediction  
5 PSIONplusm: sequential prediciton for multi-label prediction  


# Example
1. mkdir ./output/example  
2. cp example.fasta ./output/example/target.fasta  
3. python pred.py -i example -t 5  
4. result are in the ./ouput/example/output.ION (1)  output.VLG (2) output.VGS (3) output.SEQs (4) output.SEQ (5)  

# Reference  

Gao J, Cui W, Sheng Y, Ruan J, Kurgan L.PSIONplus: Accurate Sequence-Based Predictor of Ion Channels and Their Types,PLoS One, 2016, 11(4):e0152964.  

Gao J, Miao Z, Zhang Z, Wei H, Kurgan L. Prediction of Ion Channels and their Types from Protein Sequences: Comprehensive Review and Comparative Assessment. Current Drug Targets. 2019;20(5):579–592. doi:10.2174/1389450119666181022153942  

David T. Jones, Domenico Cozzetto, DISOPRED3: precise disordered region predictions with annotated protein-binding activity, Bioinformatics, Volume 31, Issue 6, 15 March 2015, Pages 857–863, https://doi.org/10.1093/bioinformatics/btu744  

Faraggi E, Zhang T, Yang Y, Kurgan L, Zhou Y. SPINE X: improving protein secondary structure prediction by multistep learning coupled with prediction of solvent accessible surface area and backbone torsion angles. Journal of computational chemistry. 2012; 33(3):259–67. Epub 2011/11/03. doi: 10.1002/jcc.21968 PMID: 22045506; PubMed Central PMCID: PMC3240697.  

Jones, D.T. (1999) Protein secondary structure prediction based on position-specific scoring matrices. J. Mol. Biol. 292:195-202.  

Altschul SF, Madden TL, Schäffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs.Nucleic Acids Res. 1997 Sep 1;25(17):3389-402.  

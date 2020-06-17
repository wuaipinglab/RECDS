		INSTALLATION AND IMPLEMENTATION OF RECDS and Third-party software

RECDS evaluate the contribution of a specific amino acid site on the HA protein in the whole history of antigenic evolution.
In RECDS, we ranked all of the HA sites by calculating the contribution scores derived from the forest of gradient boosting classifiers trained by various sequence-based and structure-based features. 

Copyright of RECDS scripts in code dir, including how to extract features, train model and color trees, are reserved by Aiping Wu lab.
When you use it you should modify python script files and set some paths based on your computing environment and third-party software installation.  

This directory contains the following subdirectories:

1. Directory/dataset/

&ensp;&ensp;&ensp;a) H3N2AntigenicCluster.strains：listing the strain names of A/H3N3 used by RECDS

&ensp;&ensp;&ensp;b) H3N2/ source: contains seven antigenic clusters derived from Smith’s studies (Koel, B.F., et al. Substitutions near the receptor binding site determine major antigenic change during influenza virus evolution. Science 2013;342(6161):976-979.) including raw HA1 protein sequences.

&ensp;&ensp;&ensp;c) H3N2/nonredundant: contains the sequences after removing redundant HA1 protein sequences and sequence alignment.

&ensp;&ensp;&ensp;d) H1N1AntigenicCluster.strains：listing the strain names of A/H1N1 used by RECDS

&ensp;&ensp;&ensp;e) H1N1/ source: contains seven antigenic clusters derived from our previous studies (Liu, M., et al. Antigenic Patterns and Evolution of the Human Influenza A (H1N1) Virus. Sci Rep 2015;5:14171.) including raw HA1 protein sequences.

&ensp;&ensp;&ensp;f) H1N1/nonredundant: contains the sequences after removing redundant HA1 protein sequences and sequence alignment.



2. Directory/code/

&ensp;&ensp;&ensp;a)Script: provides the source code of RECDS including features extracting, GBC model training and phylogenetic tree coloring.

&ensp;&ensp;&ensp;&ensp;&ensp;i.Features extracting: 


&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;1) 1_PIMAscore.py for calculating PIMA score, in which the executive command as following: 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;“./script/1_PIMAscore.py dir/dData dir/features/dPIMA”

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/dData is the virus pairs name and you can find its example format in InputH3N2Example;

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/features/dPIMA is the output feature file

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- to modify 1_PIMAscore.py file and set “seqin” path which contains virus sequences named by dData


&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;2) 2_1_run_Modeler.py and 2_2_get_pdb.py for building MODELER structure, in which the executive command 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;as following: 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;“./2_1_run_Modeler.py dir/seqName #build modeler structure;

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;./2_2_get_pdb.py dir/seqName dir/Modeler dir/newModeler #rename modeler model name;”

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- seqNmae is the virus name and you can find its example format in InputH3N2Example;

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- to modify 2_1_run_Modeler.py file and set some paths, such as runModeller, fastaDir and outDir;

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/Modeler is the output file of outDir in 2_1_run_Modeler.py;

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/newModeler is the rename output file of Modeler structure;


&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;3) 3_1_runDssp.py and 3_2_get_dssp.py for relative solvent accessibility, in which the executive command as 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;following: 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;“./3_1_runDssp.py dir/seqMap

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;./3_2_get_dssp.py dir/dData dir/features/dDssp”

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/seqMap contains the virus name and you can find its example format in InputH3N2Example;

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- to modify 3_1_runDssp.py file and set some paths which contains outdir for dssp resuls and Modeler for 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;dir/newModeler

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/features/dDssp is the output feature file

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- to modify 3_2_get_dssp.py file and set dsspfile path to outdir for dssp resuls in 3_1_runDssp.py


&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;4) 4_1_run_ddG.py and 4_2_get_ddG.py for protein stability change, in which the executive command as 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;following:

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;“./4_1_run_ddG.py dir/dData dir/dSEF

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;./4_2_get_dssp.py dir/dData dir/features/dSEF”

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/dSEF output file for ddG; 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- to modify 4_1_run_ddG.py file and set “seqin” path which contains virus sequences named by dData

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/dDssp is the output feature file

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- to modify 4_2_get_dssp.py file and set ddGDir path to outdir for SEF resuls in 4_1_run_ddG.py


&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;5) 5_1_runHSE.py and 5_2_get_HSE.py for half-sphere exposure (HSE)-up, in which the executive command as 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;following: 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;“./5_1_runHSE.py dir/dData

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;./5_2_get_HSE.py dir/dData dir/features/dHSE”

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- to modify 5_1_runHSE.py file and set some paths which contains outdir for HSE resuls and ModeDir for 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;dir/newModeler

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/features/dHSE is the output feature file

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- to modify 5_2_get_HSE.py file and set hsefile path to outdir for HSE resuls in 5_1_runHSE.py; 


&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;6) 6_1_run_Amber.py	and 6_2_get_Amber.py for Van der Waals energy and the non-polar solvation energy 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;changes, in which the executive command as following: 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;“./6_1_run_Amber.py 10000 dir/amber/10000 dir/newModeler

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;./6_2_get_Amber.py dir/dData dir/features/dAmber”

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- 10000 is modeler model name;

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/amber/10000 is output dir;

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/newModeler is modeler structure dir;

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- to modify 6_1_Amber.py file and set some paths 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- to modify 5_1_runHSE.py file and set some paths which contains outdir for HSE resuls and ModeDir for 
&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;dir/newModeler

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/features/dAmber is the output feature file;

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- to modify 6_2_get_Amber.py file and set amberfile path to outdir for amber resuls in 6_1_run_Amber.py;


&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;7) 7_merge_dfeature.py for merging six different features, in which the executive command as following:

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;“./7_merge_feature.py dir/dFeatures”

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/dFeatures is output dir;

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- to modify 7_merge_feature.py file and set featureDir path which contains above six features.


&ensp;&ensp;&ensp;&ensp;&ensp;ii. Train GBC model

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;1) 8_siteImportance.py for calculating contribution score for every amino acid position, in which the executive command as following:

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;“8_siteImportance dir/dFeatures dir/sFeatures dir/train/allImportance”

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/dFeatures is features file for antigenically different virus pairs;

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/sFeatures is features file for antigenically similar virus pairs;

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;---- dir/train/allImportance is score output file.


&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;2) 9_cross_train.py and siteImportance_diffComb.py is similar with 8_siteImportance.py and they are our 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;experiment scripts.


&ensp;&ensp;&ensp;&ensp;&ensp;iii. Custom-made phylogenetic analysis.

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;colorTree.zip：the dedicated phylogenetic visualization R package ggtree (v1.6.1) was used to plot and color 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;the tree according to the amino acid type of the candidate positions. The phylogenetic trees of A/H3N2 and 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;A/H1N1 were built under a GTR substitution model using RAxML (version 8.2.9). The GAMMA model of rate 

&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;&ensp;heterogeneity and other model parameters were estimated by RAxML.


&ensp;&ensp;&ensp;b) python: provides some scripts that will be used in the source code of RECDS (Directory/code/Script)



3. Directory/tool

&ensp;&ensp;&ensp;Directory/tool/* provides some tools that will be used in the source code of RECDS (Directory/code/Script)



Third-party software installation:  
While the majority of programs in the package 'code' are developed in the Aiping Wu lab 
here in the permission of use is released, there are some programs and databases which 
were developed by third-party groups. You may experience installation problems, please check the third-party websites.   

a) We calculated ∆∆G denoting a protein stability change upon the single point mutation through the statistical energy function (SEF).
This software is developed by Xiong Peng, and he modified application program interfaces (APIs) from original program for our use.
So, the input is very simple, which just need protein structure and mutation information (wild type, position, mutant type).
We have obtained his approval to open his modified software with our works and its software can be downloaded at https://drive.google.com/open?id=1Hls89AV6r5DNMXyDfC-qv4YhxkW_iaov.
If you use SEF, please cite as following:  

Xiong, P., et al., Protein design with a comprehensive statistical energy function and boosted by experimental selection for foldability. Nat Commun, 2014. 5: p. 5330   

b) Installing Amber  
&ensp;&ensp;fill out the form from http://ambermd.org/AmberTools14-get.html and click Download.   
&ensp;&ensp;> tar -jxvf AmberTools14.tar.bz2  
&ensp;&ensp;> rm AmberTools14.tar.bz2  
&ensp;&ensp;> cd amber14   
&ensp;&ensp;> vim ~/.bashrc  
&ensp;&ensp;copy this sentence "export AMBERHOME=$yourpath/amber14" to ~/.bashrc  
&ensp;&ensp;> source ~/.bashrc  
&ensp;&ensp;> sudo apt-get install csh flex gfortran g++ xorg-dev zlib1g-dev libbz2-dev patch python-tk python-matplotlib  
&ensp;&ensp;> ./configure gnu  
&ensp;&ensp;> source amber.sh  
&ensp;&ensp;> make install  
&ensp;&ensp;> make test  
    which will run tests and will report successes or failures

c) Installing scikit-learn  
&ensp;&ensp;i) Scikit-learn requires: Python (>= 2.6 or >= 3.3), NumPy (>= 1.6.1), SciPy (>= 0.9)
    (http://scikit-learn.org/stable/install.html);
    In the download webpage https://pypi.org/project/scikit-learn/0.18.2/, you can downlod scikit-learn-0.17.tar.gz by the link: scikit-learn-0.17.tar.gz (md5)
    or  
&ensp;&ensp;&ensp;>wget https://files.pythonhosted.org/packages/26/c2/21c612f3a1b1ba97b7b4bbd1fcdc59b475a09e25efad13fec4565ab9d563/scikit-learn-0.18.2.tar.gz   
&ensp;&ensp;ii) > tar zxvf scikit-learn-0.18.2.tar.gz  
&ensp;&ensp;&ensp;> cd scikit-learn-0.18.2  
&ensp;&ensp;&ensp;To install in your directory use:  
&ensp;&ensp;&ensp;> python setup.py install --prefix $yourpath  

d) Beyond that, you also need to install other popular softwares, such as Modeller (https://salilab.org/modeller/), dssp for SA and biopython-1.70 for import HSExposure   

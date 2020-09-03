# INSTALLATION AND IMPLEMENTATION OF RECDS and Third-party software

`RECDS` evaluate the contribution of a specific amino acid site on the HA protein in the whole history of antigenic evolution.
In RECDS, we ranked all of the HA sites by calculating the contribution scores derived from the forest of gradient boosting classifiers trained by various sequence-based and structure-based features. 

Copyright of `RECDS` scripts in code dir, including how to extract features, train model and color trees, are reserved by Aiping Wu lab.

***Important: when you use it you should modify python script files and set some paths based on your computing environment and third-party software installation.***

The repo contains the following directories:


## 1. Raw data used in this project

The antigenic clustering and raw sequence used in `RECDS` can be found in `dataset/` directory.

### H3N2

`H3N2AntigenicCluster.strains`: The strain names of `A/H3N2` used by `RECDS`

`H3N2/source/`: Seven antigenic clusters derived from Smith’s studies including raw HA1 protein sequences.

`H3N2/nonredundant/`: The sequences after removing redundant HA1 protein sequences and sequence alignment.

_Koel, B.F., et al. Substitutions near the receptor binding site determine major antigenic change during influenza virus evolution. Science 2013;342(6161):976-979._

### H1N1

`H1N1AntigenicCluster.strains`: The strain names of `A/H1N1` used by `RECDS`

`H1N1/source/`: Seven antigenic clusters derived from our previous studies including raw HA1 protein sequences.

`H1N1/nonredundant/`: The sequences after removing redundant HA1 protein sequences and sequence alignment.

_Liu, M., et al. Antigenic Patterns and Evolution of the Human Influenza A (H1N1) Virus. Sci Rep 2015;5:14171._


## 2. Pipeline and result of the project

Source code of the `RECDS` pipeline can be found in `code/sourceCode/script/`. The major steps includes a) features extracting, b) GBC model training and c) phylogenetic tree coloring. The dependency scripts can be found in `code/sourceCode/python/` and `tools/`.

### a) Features extracting:

Including the PIMA score and other metrics calculated from predicted 3D structure.

#### PIMA score

1. Modify `1_PIMAscore.py` file and set “seqin” path which contains virus sequences named by `dData`
2. `dir/dData` is the virus pairs name and you can find its example format in `InputH3N2Example`
3. `dir/features/dPIMA` is the output feature file

```bash
./1_PIMAscore.py dir/dData dir/features/dPIMA
```

#### MODELER structure

1. Modify `2_1_run_Modeler.py` file and set some paths, such as runModeller, fastaDir and outDir
2. `seqName` is the virus name and you can find its example format in `InputH3N2Example`
3. `dir/Modeler` is the output file of outDir in `2_1_run_Modeler.py`
4. `dir/newModeler` is the rename output file of Modeler structure

```bash
./2_1_run_Modeler.py dir/seqName # build modeler structure
./2_2_get_pdb.py dir/seqName dir/Modeler dir/newModeler # rename modeler model name
```

#### Relative solvent accessibility

1. Modify `3_1_runDssp.py` file and set some paths which contains outdir for dssp resuls and Modeler for `dir/newModeler`; 
2. `dir/seqMap` contains the virus name and you can find its example format in `InputH3N2Example`
3. Modify `3_2_get_dssp.py` file and set dsspfile path to outdir for dssp resuls in `3_1_runDssp.py`
4. `dir/features/dDssp` is the output feature file

```bash
./3_1_runDssp.py dir/seqMap
./3_2_get_dssp.py dir/dData dir/features/dDssp
```

#### Protein stability change

1. Modify `4_1_run_ddG.py` file and set “seqin” path which contains virus sequences named by `dData`
2. `dir/dSEF` output file for ddG
3. Modify `4_2_get_ddG.py` file and set ddGDir path to outdir for SEF resuls in `4_1_run_ddG.py`
4. `dir/dSEF` is the output feature file

```bash
./4_1_run_ddG.py dir/dData dir/dSEF
./4_2_get_ddG.py dir/dData dir/features/dSEF
```

#### Half-sphere exposure (HSE)-up

1. Modify `5_1_runHSE.py` file and set some paths which contains outdir for HSE resuls and ModelDir for `dir/newModeler`
2. Modify `5_2_get_HSE.py` file and set hsefile path to outdir for HSE resuls in `5_1_runHSE.py`
3. `dir/features/dHSE` is the output feature file

```bash
./5_1_runHSE.py dir/dData
./5_2_get_HSE.py dir/dData dir/features/dHSE
```

#### Energy changes

1. Modify `6_1_run_Amber.py` file and set some paths
2. `10000` is modeler model name
3. `dir/amber/10000` is output dir
4. `dir/newModeler` is modeler structure dir
5. Modify `6_2_get_Amber.py` file and set amberfile path to outdir for amber resuls in `6_1_run_Amber.py`
6. `dir/features/dAmber` is the output feature file

```bash
# Van der Waals energy and the non-polar solvation energy
./6_1_run_Amber.py 10000 dir/amber/10000 dir/newModeler
./6_2_get_Amber.py dir/dData dir/features/dAmber
```

#### Merge six different features

1. Modify `7_merge_feature.py` file and set featureDir path which contains above six features
2. `dir/dFeatures` is output dir

```bash
./7_merge_feature.py dir/dFeatures
```

### b) Train GBC model

Use the features extracted in the previous steps to train Gradient Boost Classifier model.

1. `9_cross_train.py` and `siteImportance_diffComb.py` is similar to `8_siteImportance.py` and they are our experiment scripts
2. `dir/dFeatures` is features file for antigenically different virus pairs
3. `dir/sFeatures` is features file for antigenically similar virus pairs
4. `dir/train/allImportance` is score output file

```bash
 # calculate contribution score for every amino acid position
./8_siteImportance dir/dFeatures dir/sFeatures dir/train/allImportance
```

### c) Custom-made phylogenetic analysis.

Build phylogenetic tree, plot and color the tree according to the amino acid type of the candidate positions

1. The phylogenetic trees of `A/H3N2` and `A/H1N1` were built under a GTR substitution model using `RAxML` (version 8.2.9). The GAMMA model of rate heterogeneity and other model parameters were estimated by `RAxML`.
2. The dedicated phylogenetic visualization R package `ggtree` (v1.6.1) was used to make the plot.
3. The result can be found in `colorTree.zip`.


## 3. Third-party software installation:  
While the majority of programs in the package `code/sourceCode/` are developed in the Aiping Wu lab here in the permission of use is released, there are some programs and databases which were developed by third-party groups. You may experience nstallation problems, please check the third-party websites.   

### Install SEF
We calculated ∆∆G denoting a protein stability change upon the single point mutation through the statistical energy function (SEF). This software is developed by Xiong Peng, and he modified application program interfaces (APIs) from original program for our use. So, the input is very simple, which just need protein structure and mutation information (wild type, position, mutant type). We have obtained his approval to open his modified software with our works and its software can be downloaded at https://drive.google.com/open?id=1Hls89AV6r5DNMXyDfC-qv4YhxkW_iaov. If you use SEF, please kindly cite as following:  

_Xiong, P., et al., Protein design with a comprehensive statistical energy function and boosted by experimental selection for foldability. Nat Commun, 2014. 5: p. 5330_

### Install Amber  
Fill out the form from http://ambermd.org/GetAmber.php and click Download.   
```bash
tar -jxvf AmberTools14.tar.bz2
rm AmberTools14.tar.bz2
cd amber14
echo "export AMBERHOME=$yourpath/amber14" >> ~/.bashrc
source ~/.bashrc
sudo apt-get install csh flex gfortran g++ xorg-dev zlib1g-dev libbz2-dev patch python-tk python-matplotlib
./configure gnu
source amber.sh
make install
make test  # which will run tests and will report successes or failures
```

### Install scikit-learn  
`scikit-learn` requires: `Python` (>= 2.6 or >= 3.3), `NumPy` (>= 1.6.1), `SciPy` (>= 0.9); Here is one way to install
```bash
wget https://files.pythonhosted.org/packages/26/c2/21c612f3a1b1ba97b7b4bbd1fcdc59b475a09e25efad13fec4565ab9d563/scikit-learn-0.18.2.tar.gz
tar zxvf scikit-learn-0.18.2.tar.gz
cd scikit-learn-0.18.2
python setup.py install --prefix $yourpath # To install in your directory use
```
### Other dependencies
Beyond that, you also need to install other popular softwares, such as `Modeller` (https://salilab.org/modeller/), `dssp` for SA and `biopython` for `HSExposure` module.

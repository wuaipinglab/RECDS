		INSTALLATION AND IMPLEMENTATION OF RECDS and Third-party software

Copyright of RECDS scripts in code dir,including how to extract features, train model and color trees, are reseved by Aiping Wu lab.
When you use it you should modify python script files and set some paths based on your computing environment and third-party software installation.  


Third-party software installation:  
While the majority of programs in the package 'code' are developed in the Aiping Lab 
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
&ensp;&ensp;>wget https://files.pythonhosted.org/packages/26/c2/21c612f3a1b1ba97b7b4bbd1fcdc59b475a09e25efad13fec4565ab9d563/scikit-learn-0.18.2.tar.gz   
&ensp;&ensp;ii) > tar zxvf scikit-learn-0.18.2.tar.gz  
&ensp;&ensp;&ensp;&ensp;> cd scikit-learn-0.18.2  
&ensp;&ensp;&ensp;&ensp;To install in your directory use:  
&ensp;&ensp;&ensp;&ensp;> python setup.py install --prefix $yourpath  

d) Beyongd that, you also need to install other popular softwars, such as Modeller (https://salilab.org/modeller/), dssp for SA and biopython-1.70 for import HSExposure   

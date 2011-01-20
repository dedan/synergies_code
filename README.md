my work in yifat's lab
======================

the first commit to this repository is the the last commit
of the svn repository that I used during my BA thesis

External Toolboxes which are required to run this project:

* [circular statistics](http://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics)
* [ica from eeglab](http://sccn.ucsd.edu/eeglab)


How to use these scripts
========================

Data collection and preprocessing is done by the files in the preprocessing folder. 
There are two *callable* files. 

* **runscript4emgdat** to create the emgdata (in folder EMGdat)
* **runscript4evoked** to compute the average windows.

The data computed in this first steps can then be used in the analysis files 

* **compute_synergies** to compute the synergies and do residual tests on the data from natural movement
* --> this will result in all_data_monkname files to be used in nonevoked_analysis


Explanation of folder content and how to use the files:
-------------------------------------------------------

* config

    contains all the configuration files
    
* data_validation

    these files are used to assure the correctness of data or to examine raw data
    
* method_validation

    investigate the quality and suitability of methods, e.g. how stable is the 
    nmf factorization, or: does pcaica come to the same results
    
* lib

    helper methods and functions
    
* locations

    investigate whether the synergies have a topographical organization in the cortex
    
* nonevoked_syns

    Analysis whether we find muscle synergies during natural movement of the monkey
    
* preprocess

    Preprocessing of the raw data. Read emg files, compute the responses or average
    windows, and so on
    
* synergists

    analysis of synergies of responses evoked by cortical stimulation
    
* test_data 

    data which is used in the data and method validation functions

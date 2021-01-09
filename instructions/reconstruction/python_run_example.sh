
############ INSTALL DEPENDENCIES ############


# Step 1: install conda (miniconda)
# E.g. 
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
# you can view more instructions here 
# https://docs.conda.io/en/latest/miniconda.html#linux-installers

# Step 2: create new conda environment with python 2 dependency 
conda create --name py38 python=3.8
conda activate py38

# Step 3: install all python dependencies for reconstruction pipeline via pip 
cd reconstruction_python/
pip install -r /requirements.txt
# see footnote about `pip`



############ RUN EXAMPLE ############

# Step 1: 
# open reconstruction_python/example.sh and substitute the correct `MRN` and `name` from `input/subjectList_MRU.csv`. Remember to never store PHI data on github.

# Step 2: 
bash example.sh

# NOTES 
# it is recommended that you run reconstruction on a server with lots of memory... e.g. no less than 100Gb RAM 
# so instead of running it on your local machine you can ssh into one of the following CRL machines and run the reconstruction. e.g 
ssh boreas #or ssh ganymede/auster/zephyr
cd reconstruction_python/
bash example.sh



############ REQUIRED CHANGES FOR PROCESSING NEW RECONSTRUCTIONS ############

# Step 1: 
# To launch a new reconstruction, make the appropriate changes in example.sh


# NOTES 
# 1. Please note that the specified inputs args to `processSubjectRLTI.py` are the default recon parameters and had been proven to work well. 
# 2. We run the example reconstruction with `no-bmd` flag, which turns off bulk motion detection. This makes the python recon most similar to the matlab reconstruction. However, removing this flag may improve the reconstruction result. 
# 3. Further note that the reconstruction can be run with more sophisticated motion correction. The set of libraries specified in `reconstruction_python/ext` can be used to this extent. 
# 4. I had REMOVED a lot of files from `reconstruction_python/ext` as they had contained PHI information. These files are not necessary to start reconstruction, however, if you want to get hold of them - please refer to the notes left in each file - which will specify the location of the original libraries on the CRL filesystem. 











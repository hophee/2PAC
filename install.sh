#!/usr/bin/bash
eval "$(conda shell.bash hook)"

conda env create -f environment.yml -y
#PRIMER3
git clone https://github.com/primer3-org/primer3.git primer3
cd primer3/src
make
cd ../../

#virtualPCR
git clone https://github.com/rkalendar/virtualPCR.git

#CHOPCHOP
git clone https://bitbucket.org/valenlab/chopchop.git
#change config pathes

#callPrimer3 for R
wget https://gist.githubusercontent.com/IdoBar/5e78ae7a5cc7277a04b126ce6f595d6e/raw/45c60662f3479f41765bce839835c4988a7e5b36/callPrimer3.R


# rlProtGen (rL protein generator)
Generation of BackBones conditioned on theozyme geometry

Owners: Yakov Kipnis and Ariel Ben-Sasson

Version 0: unsupervised RL for generation of protein backbone from a set of constrained rotamers constituting a theozyme geometry

The goal here is to explore idea that protein backbones conditioned on external structural constraints represented by theozymes can be generated using trained RL agent from limited set of structural "letters".

Letters are fragments of native protein backbones, representing recurrent structural motifs.
There are two kinds of letters:
1. 4mer
2. 6mer

4mers originate from continuous fragments of native protein backbones. Fragments need to be as short as possible to minimize their diversity and allow their clustering into finite number of structural "letters", while allowing unambiguous assembly into larger structures by simple geometrical superposition of letters 

Instructions to build the conda env:
1.  clone environment:
    on a linux system:  conda env create -f environment.yml
    on a non linux system: conda env create -f environmentMin.yml
2.  activate environment:
    conda activate rlProtGen
3.  stage environemnt kernel to use on a jupyter notebook:
    python -m ipykernel install --user --name rlPrtoGen --display-name=rlPrtoGen
4.  check your kernels:
    jupyter kernelspec list
5.  open the notebook, select the "rlProtGen" kernel and enjoy. (or not)

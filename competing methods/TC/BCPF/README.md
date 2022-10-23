# Bayesian-CP-Factorization (BCPF)
Matlab code of Bayesian CP Factorization for Tensor Completion


We provide two demo codes:
I.  DemoBayesCP.m:        Demonstration on synthesic data
II. DemoBayesCP_Image.m   Demonstration for image completion


The package includes four algorithms:
1. BCPF.m           BCPF for fully observed tensor
2. BCPF_TC.m        BCPF for incomplete tensor 
3. BCPF_IC.m        BCPF for image completion
4. BCPF_MP.m        BCPF using mixture priors for image completion


In this package, we used the tensor toolbox 2.5, which is downloaded from (http://www.sandia.gov/~tgkolda/TensorToolbox)

The tools for visualization of tensor with voxels is from Tensorlab (http://www.tensorlab.net/)

Referecne:
  Q. Zhao, L. Zhang, and A. Cichocki. Bayesian CP factorization of incomplete tensors with automatic rank determination. IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI), 37(9):1751–1763, 2015.

Outline:

    Folder structure
    Usage
    Citation
    URL of third-party toolboxes and functions


Folder structure:

    data\           : example data 
        |-- testMSI_1.mat      : example MSI data 1
        |-- testMSI_2.mat      : example MSI data 2
        |-- testVideo.mat      : example vedio data 
        |-- groundtruth_at_16th_frame.bmp   : groundtruth of vedio data at 16th frame
    lib\            : library of code

        |-- tensor_toolbox\    : toolbox for tensor operations [1]
        |-- my_tensor_toolbox\ : another toolbox for tensor operations 
        |-- KBR\               : subfunctions of KBR method
        |-- quality_assess\    : functions of quality assessment indices
        |       |-- ssim_index.m            : SSIM [2]
        |       |-- FeatureSIM.m            : FSIM [3]
        |       |-- ErrRelGlobAdimSyn.m     : ERGAS
        |       |-- SpectAngMapper.m        : SAM
        |       |-- MSIQA.m                 : interface for calculating PSNR and the four indices above
        |-- compete_methods : competing methods (down loaded or implemented based on reference papers)
        |       |-- HOrpca\                 : HOrpca toolbox [4]
        |       |-- inexact_alm_rpca\       : matrix ADMM-based rpca toolbox [5]
        |       |-- tensor_SVD\             : tSVD toolbox [6]
        |       |-- TMac\                   : TMac toolbox [7]
        |       |-- TenCompletion\          : McpTC and ScadTC toolbox [8]
        |       |-- LRTV_TC.m               : functions for Trace/TV method [9]
        |       |-- Mc_adm.m                : functions for MC_ALM method [10]
        |       |-- Tc_liu.m                : functions for HaLRTC method [11]
        |-- togetGif.m         : function for getting GIF result
        |-- showVideoResult.m  : scripts that show the result of background subtraction experiment
        |-- showMSIResult.m    : scripts that show the result of MSI completion experiment
    Demo_TC_MSI.m     : scripts that applies the tensor-completion methods and calculates the QA indices
    Demo_RPCA_video.m : scripts that applies the robust-principal-component-analysis methods and calculates the QA indices
    KBR_TC.m          : function for sloving the intrinsic tensor sparsity based tensor completion model 
    KBR_RPCA.m        : core function of the the intrinsic tensor sparsity based tensor RPCA model 
 

Usage:
    
  (1) For MSI completion experiment, you can simply follow these steps:
        1.Re-arrange the MSI into [0, 1].
        2.Add the folder 'lib'into path, and use the function KBR_TC as follows:
            [ restored_img ] = KBR_TC( courrupted_img, Omega)
    Please type 'help KBR_TC ' to get more information.

    You may find example codes in file Demo_TC_MSI.m

    Also, you can use the demo to see some comparison. You can:
      1. Type 'Demo_TC_MSI' to to run various methods and see the pre-computed results.
      2. Use 'help Demo_TC_MSI' for more information.
      3. Change test MSI by simply modifying variable 'dataname' in Demo_TC_MSI.m (NOTE: make sure your MSI
         meets the format requirements).
      4. Change sampling rate by modifying variables  'sample_ratio ' in Demo_TC_MSI.m
      5. Select competing methods by turn on/off the enable-bits in Demo_TC_MSI.m

  (2) For background subtraction experiment, you can simply follow these steps:
        1.Re-arrange the video into [0, 1].
        2.Use the function KBR_RPCA as follows:
            [ restored_img ] = KBR_RPCA( tensorData,beta,gamma)
    Please type 'help KBR_RPCA ' to get more information.

    You may find example codes in file Demo_RPCA_video.m

    Also, you can use the demo to see some comparison. You can:
      1. Type 'Demo_RPCA_video' to to run various methods and see the pre-computed results.
      2. Use 'help Demo_RPCA_video' for more information.
      3. Change test video data by simply modifying variable 'dataname' in Demo_RPCA_video.m (NOTE: make sure your video
         meets the format requirements).
      4. Select competing methods by turn on/off the enable-bits in Demo_RPCA_video.m

Citation:

    Qi Xie, Qian Zhao, Deyu Meng*, & Zongben Xu
    Kronecker-Basis-Representation Based Tensor Sparsity and Its Applications to Tensor Recovery[J]. 
    IEEE Transactions on Pattern Analysis & Machine Intelligence, 2017, PP(99):1-1. (accepted)

    BibTeX:
      @article{Qi2017Kronecker,  
      title={Kronecker-Basis-Representation Based Tensor Sparsity and Its Applications to Tensor Recovery},  
      author={Qi, Xie and Qian, Zhao and Meng, Deyu and Xu, Zongben},  
      journal={IEEE Transactions on Pattern Analysis & Machine Intelligence},  
      volume={PP},  
      number={99},  
      pages={1-1},  
      year={2017},
      }

URL  of the toolboxes and functions and citation of competing methods:

    [1]  tensor_toolbox     http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.5.html
    [2]  ssim_index.m       https://ece.uwaterloo.ca/~z70wang/research/ssim/
    [3]  FeatureSIM.m       http://www4.comp.polyu.edu.hk/~cslzhang/IQA/FSIM/FSIM.htm
    [4]  HOrpca             https://sites.google.com/site/tonyqin/research
    [5]  inexact_alm_rpca   http://perception.csl.illinois.edu/matrix-rank/home.html
    [6]  t-SVD              http://www.ece.tufts.edu/~shuchin/software.html
    [7]  Tmac               http://www.caam.rice.edu/~yx9/TMac/
    [8]  W. Cao, Y. Wang, C. Yang, X. Chang, Z. Han, and Z. Xu. Foldedconcave penalization approaches to tensor completion.
         Neurocomputing, 152:261¨C273, 2015.
    [9]  M. Golbabaee and P. Vandergheynst. Joint trace/tv norm minimization: A new efficient approach for spectral compressive
         imaging. pages 933¨C936, 2012.
    [10] Z. Lin, M. Chen, and Y. Ma. The augmented lagrange multiplier method for exact recovery of corrupted low-rank matrices. 
         Eprint Arxiv, 9, 2010.
    [11] J. Liu, P. Musialski, P. Wonka, and J. Ye. Tensor completion for estimating missing values in visual data. IEEE Trans. 
         Pattern Analysis & Machine Intelligence, 35(1):208¨C220, 2013.
    




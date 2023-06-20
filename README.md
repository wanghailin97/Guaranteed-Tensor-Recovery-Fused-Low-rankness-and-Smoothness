# Guaranteed-Tensor-Recovery-Fused-Low-rankness-and-Smoothness
The released code of t-CTV algorithms for tensor completion (TC) and tensor robust principal componnet analysis (TRPCA), mainly proposed in the paper "Guaranteed Tensor Recovery Fused Low-rankness and Smoothness" which is accepted by the top journal T-PAMI, 2023.

Citation: 

@ARTICLE{10078018,

  author={Wang, Hailin and Peng, Jiangjun and Qin, Wenjin and Wang, Jianjun and Meng, Deyu},
  
  journal={IEEE Transactions on Pattern Analysis and Machine Intelligence}, 
  
  title={Guaranteed Tensor Recovery Fused Low-rankness and Smoothness}, 
  
  year={2023},
  
  pages={1-17},
  
  doi={10.1109/TPAMI.2023.3259640}}


------ Description-----

--1. The proposed t-CTV regularizer is very simple and powerful, which can simultaneously encode the global low-rankness and local smoothness of a tensor data by using a single regularizer. This is totally different with current joint low-rank and smooth modeling manner.

--2. The t-CTV induced TC and TRPCA problems are proven to be exactly recovered. As far as we known, this should be the first theoretical gaurantee for joint low-rank and smooth tensor recovery. 

--3. The very impressive result is that the t-CTV based TC model has a lower bound of sampling complexity than purely low-rank regularizer based models, purely smooth regularizer based models, and low-rank regularizer plus smooth regularizer based models. This explains that the proposed t-CTV based TC method can still estimate partial information in tensor completion tasks under very high missing rates, such as 99%, 99.5%. 

--4. The t-CTV model is totally parameter-free and beats a lot of classic and efficient tensor recovery methods, including SNN(HaLRTC), TNN, BCPF, KBR, SPC-TV, SNN-TV, TNN-TV and so on. The released code gives a very convenient demo test. 


------May this work is helpful to you!-----

-------------------------------------------wrriten by Hailin Wang, 2023-6-20--------------------------------------

--if you have any question, please contact me without hesitation--email: wanghailin97@163.com





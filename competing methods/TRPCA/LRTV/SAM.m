function SAM_value=SAM(T,H)
SigmaTR=T*H';
SigmaT2=T*T';
SigmaR2=H*H';
SAM_value=acosd(SigmaTR/sqrt(SigmaT2*SigmaR2));
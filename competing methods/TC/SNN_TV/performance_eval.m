function [perform_RSE,perform_PSNR] = performance_eval()
    filename=cell(8,1);
    filename{1}='Result/airplane';
    filename{2}='Result/baboon';
    filename{3}='Result/barbara';
    filename{4}='Result/facade';
    filename{5}='Result/house';
    filename{6}='Result/lena';
    filename{7}='Result/peppers';
    filename{8}='Result/sailboat';
    
    perform_RSE=zeros(2,8);
    perform_PSNR=zeros(2,8);
  %  for i=1:7 % seven methods
        for k=1:8 % missing rate
            %for j=1:8 % images
             for j=1:8
                myname=[filename{j} '_Data.mat'];
                load(myname);
                myname=[filename{j} '_Result.mat'];
                load(myname);

                
                
                tmpRSE_LRTC_TV1=rse(A,myResult{1,k});
                tmpPSNR_LRTC_TV1=psnr(A,myResult{1,k});
                perform_RSE(1,k)=perform_RSE(1,k) + tmpRSE_LRTC_TV1;
                perform_PSNR(1,k)=perform_PSNR(1,k) + tmpPSNR_LRTC_TV1;
                
                tmpRSE_LRTC_TV2=rse(A,myResult{2,k});
                tmpPSNR_LRTC_TV2=psnr(A,myResult{2,k});
                perform_RSE(2,k)=perform_RSE(2,k) + tmpRSE_LRTC_TV2;
                perform_PSNR(2,k)=perform_PSNR(2,k) + tmpPSNR_LRTC_TV2;
            end
        end
         perform_RSE=perform_RSE/8;
         perform_PSNR=perform_PSNR/8;
         figure(1);plot(perform_RSE');
         figure(2);plot(perform_PSNR');
        
end


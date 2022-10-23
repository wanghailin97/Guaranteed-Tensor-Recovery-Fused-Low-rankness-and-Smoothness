function [ ] = test_performance()
    filename=cell(8,1);
    filename{1}='airplane';
    filename{2}='baboon';
    filename{3}='barbara';
    filename{4}='facade';
    filename{5}='house';
    filename{6}='lena';
    filename{7}='peppers';
    filename{8}='sailboat';
    
    for ff=1:8
        myName=sprintf('TestImages/%s.bmp',filename{ff});
        A=imread(myName);
        figure(1);  imshow(A);

        myrate=0.60:0.05:0.95;
        myResult=cell(2,numel(myrate));
        A=double(A)/255.0;

        for iterate=1:numel(myrate)
            rate=1 - myrate(iterate);
            [row, col, channel]=size(A);
            B=zeros([row, col, channel]);
            mark=true([row, col, channel]);

            counter=1;
            for i=1:row
                for j=1:col
                    for k=1:channel
                        if(rand()<rate)
                            index(counter,1)=i;
                            index(counter,2)=j;
                            index(counter,3)=k;
                            value(counter)=A(i,j,k);
                            B(i,j,k)=A(i,j,k);
                            mark(i,j,k)=false;
                            counter=counter+1;
                        end
                    end
                end
            end
            figure(2);imshow(B);

            tsize=[row, col, channel];
            N=3;
            lambda=0.02;
            alpha=[1/N, 1/N, 1/N];
            beta=[1,1,0];

            fprintf('------TR_LRTV---------- \n');
            Z_TRLRTV=LRTC_TV_I(index, value, lambda, alpha, beta, tsize, N );
            figure(3);imshow(Z_TRLRTV);

            fprintf('--------------TR_LRTV2-------------------\n');
            lambda_1=0.5;
            lambda_2=1000;

            Z_TRLRTV2=LRTC_TV_II(index, value, lambda_1, lambda_2 ,alpha, beta, tsize, N );

            myResult{1,iterate}=Z_TRLRTV;  %LRTC_TV_I
            myResult{2,iterate}=Z_TRLRTV2; %LRTC_TV_II

            myName=['Result/' filename{ff} '_Result.mat'];
            save(myName,'myResult','myrate');

            myName=['Result/' filename{ff}, '_Data.mat'];
            save(myName,'A');
        end
    end
    
    [perform_RSE,perform_PSNR] = performance_eval();
    figure(1);plot(perform_RSE);
    figure(2);plot(perform_PSNR);
end


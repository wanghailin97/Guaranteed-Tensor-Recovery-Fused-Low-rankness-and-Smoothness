function showVideoResult(D,B, methodname,enList,nfram)
%% Use to show the video result of RPCA methods

numLine = length(enList);
% figure('units','normalized','position',[0.05,0.482-0.29*3/2,0.9,0.29*3],'name','video');
%% UI
close all;
figure('units','normalized','position',[0.05,0.482-0.29*3/2,0.9,0.29*3],'name','video');
sld = uicontrol('Style', 'slider',...
    'Min',1,'Max',50,'Value',nfram,...
    'Position', [220 20 500 20],...
    'Callback', @TheFram);
showShow(nfram)

    function showShow(nfram)
        txt = uicontrol('Style','text',...
            'Position',[220 45 500 20],...
            'String',['Showing the ', num2str(nfram) ,'th fram']);
        set(txt,'Fontsize',13)
        tempD = D(:,:,nfram);
        minI = min(tempD(:));
        maxI = max(tempD(:));
        ImD      = (tempD-minI)/(maxI-minI);
        for j = 1:numLine
            clear Par
            Par.Osz    = [0.15,0.15];
            Par.lineW  = 2;
            Par.times  = 3;
            Par.outX   = 0.0;
            Par.ifDB   = 1;
            %         Par.maxP   = 0.162;
            %         Par.minP   = 0.055;
            Par.color  = [1,0.2,0.2];
            tempB = B{enList(j)}(:,:,nfram);
            ImB      = (tempB-minI)/(maxI-minI);
            Y = WindowBig(ImB,[0.48,0.19],Par); % WindowGig(X,local,Osz,times,outX)
            %     sizeH = size(Y);
            subplot(3, numLine ,j); imshow(ImD); title(' Original Video')
            subplot(3,numLine ,j+numLine);
            imshow(Y);  title( [methodname{enList(j)}, ' background']);
            
            Par.Osz    = [0.15,0.15];
            Par.color  = [0.2,1,0.2];
            Par.ifDB   = 1;
            Par.maxP   = 0.4;
            Par.minP   = 0.0;
            
            ImF  = abs(ImD - ImB);
            Y = WindowBig(ImF,[0.58,0.3],Par); % WindowGig(X,local,Osz,times,outX)
            subplot(3,numLine ,j+numLine*2);
            imshow(Y); title( [methodname{enList(j)}, ' foreground']);
        end
    end

    function TheFram(source,callbackdata)
        nfram =  round(source.Value);
        showShow(nfram)
    end
end
function showMSIResult(Re_msi,T,T_miss, methodname,enList,band2show)
%% Use to show the video result of TC methods

numLine = ceil((length(enList)+2)/5);
% figure('units','normalized','position',[0.05,0.482-0.29*3/2,0.9,0.29*3],'name','video');
%% UI
close all;
figure('units','normalized','position',[0.05,0.482-0.29*numLine/2,0.9,0.29*numLine],'name','MSI');
sld = uicontrol('Style', 'slider',...
    'Min',1,'Max',31,'Value',band2show,...
    'Position', [220 20 500 20],...
    'Callback', @TheFram);
showShow(band2show)

    function showShow(band2show)
        txt = uicontrol('Style','text',...
            'Position',[220 45 500 20],...
            'String',['Showing the ', num2str(band2show) ,'th fram']);
        set(txt,'Fontsize',13)
        tempT = T(:,:,band2show);
        minI = min(tempT(:));
        maxI = max(tempT(:));
        Par.Osz    = [0.15,0.15];
        Par.lineW  = 2;
        Par.times  = 2.5;
        Par.outX   = 0.0;
        Par.ifDB   = 1;
        local      = [0.5,0.35];
        %         Par.maxP   = 0.162;
        %         Par.minP   = 0.055;
        Par.color  = [1,0.2,0.2];
        ImB      = (T(:,:,band2show)-minI)/(maxI-minI);
        Y = WindowBig(ImB,local,Par);
        subplot(numLine,5,1); imshow(Y),title( 'Original');
        ImB      = (T_miss(:,:,band2show)-minI)/(maxI-minI);
        Y = WindowBig(ImB,local,Par);
        subplot(numLine,5,2); imshow(Y),title( 'Corrupted');
        for i = 1:length(enList)
            ImB      = (Re_msi{enList(i)}(:,:,band2show)-minI)/(maxI-minI);
            Y = WindowBig(ImB,local,Par);
            subplot(numLine,5,i+2);
            imshow(Y);title( methodname{enList(i)});
        end
    end


    function TheFram(source,callbackdata)
        band2show =  round(source.Value);
        showShow(band2show)
    end
end
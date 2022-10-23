function showHSI(hsi,band2show)
%% Use to show an HSI 
close all;

if nargin < 2
    band2show = 7;
end
bandnum = size(hsi, 3);

figure('units','normalized','position',[0.05,0.482-0.29/2,0.9,0.29], 'name','HSI');
sld = uicontrol('Style', 'slider',...
    'Min',1,'Max',bandnum,'Value',band2show,...
    'Position', [645 20 500 20],...
    'Callback', @TheFram);
showShow(band2show)

    function showShow(band2show)
        txt = uicontrol('Style','text',...
            'Position',[645 45 500 20],...
            'String',['Showing the ', num2str(band2show) ,'th band']);
        set(txt,'Fontsize',13)
        imshow(hsi(:,:,band2show)); 
    end

    function TheFram(source,callbackdata)
        band2show =  round(source.Value);
        showShow(band2show)
    end
end
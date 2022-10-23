function [Y, maxP, minP] = WindowBig(X,local,Par)

if nargin<2
    local=[0.1,0.1];
end
if nargin<3
    Osz    = [0.1,0.1];
    outX   = 0.1;
    lineW  = 1;
    times  = 2;
    color  = [1,0,0];
    ifDB   = 0;
else
    if isfield(Par,'Osz')
        Osz = Par.Osz;
    else
        Osz = [0.1,0.1];
    end
    if isfield(Par,'outX')
        outX = Par.outX;
    else
        outX = 0.1;
    end
    if isfield(Par,'lineW')
        lineW = Par.lineW;
    else
        lineW = 1;
    end
    if isfield(Par,'times')
        times = Par.times;
    else
        times = 2;
    end
    if isfield(Par,'color')
        color = Par.color;
    else
        color = [1,0,0];
    end
    if isfield(Par,'ifDB')
        ifDB = Par.ifDB;
    else
        ifDB = 0;
    end
    
end



sizeX   = size(X);
sizeP   = round(sizeX(1:2).*Osz);
sizeOut = round(sizeX(1:2)*outX);
sizeY   = [sizeX(1:2)+2*sizeOut,3];

x1      = max(round(sizeX(1:2).*local),1);
x2      = min(x1+sizeP,sizeX(1:2));

Y       = ones(sizeY);
if length(sizeX)==2
    for i = 1:3
        Y(sizeOut(1)+1:sizeOut(1)+sizeX(1),sizeOut(2)+1:sizeOut(2)+sizeX(2),i) = X;
    end
elseif length(sizeX)==3
    Y(sizeOut(1)+1:sizeOut(1)+sizeX(1),sizeOut(2)+1:sizeOut(2)+sizeX(2),1:3) = X;
else
    error('Please input a image data')
end
Patch   = Y(x1(1):x2(1),x1(2):x2(2),1:3);
if ifDB
    if isfield(Par,'maxP');maxP=Par.maxP;else maxP  = max(Patch(:));end
    if isfield(Par,'minP');minP=Par.minP;else minP  = min(Patch(:));end
    Patch = (Patch-minP)/(maxP-minP); 
end
rePatch = imresize(Patch,times);

resizeP = size(rePatch);
y2      = sizeY(1:2);
y1      = max(1,y2-resizeP(1:2)+1);
Y(y1(1):y2(1),y1(2):y2(2),1:3)=rePatch;

for i = 1:3
        Y(x1(1):x2(1),[x1(2):x1(2)+lineW-1,x2(2)-lineW+1:x2(2)],i)=color(i);
        Y([x1(1):x1(1)+lineW-1,x2(1)-lineW+1:x2(1)],x1(2):x2(2),i)=color(i);

        Y(y1(1):y2(1),[y1(2):y1(2)+lineW-1,y2(2)-lineW+1:y2(2)],i)=color(i);
        Y([y1(1):y1(1)+lineW-1,y2(1)-lineW+1:y2(1)],y1(2):y2(2),i)=color(i);

end

for i = [2,3]
    
end
end
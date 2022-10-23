function [SigmaX,svp]=ClosedWNNM(SigmaY,C,oureps)
% solving the following problem
%         sum(log(SigmaY+oureps))+1/2*||Y-X||_F^2
% where oureps is a small constant

temp   = (SigmaY-oureps).^2-4*(C-oureps*SigmaY);
ind    = find (temp>0);
svp    = length(ind);
SigmaX = max(SigmaY(ind)-oureps+sqrt(temp(ind)),0)/2;
end
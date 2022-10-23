
function [Core, U] = tensorSVD(D,rank)
% performing Higth order SVD
dim        = length(size(D));
sizeD      = size(D);
U          = cell(dim,1);
Core       = D;
if nargin<2
    rank = sizeD;
end

for i = 1:dim
    UnfoD = Unfold( double(D), sizeD, i);
    if size(UnfoD,1)<size(UnfoD,2)
        DxD   = UnfoD*UnfoD';  
        if sizeD(i) <= rank(i)
            [temp , ~, ~]    = svd(DxD);
        else
            eigsopts.disp = 0;
            [temp , ~]    =  eigs(DxD, rank(i) , 'LM', eigsopts);
        end
        U{i}          = temp;
        Core          = ttm(tensor(Core), temp', i);
    else
        [temp , ~, ~]    = svd(UnfoD,'econ');
        U{i}          = temp(:,1:rank(i));
        Core          = ttm(tensor(Core), temp', i);  
    end
end

Core = double(Core);
end


function [ outpatch ] = imageTopatch( image,par )
% this function is to slide the 3-D image into patches
% -------input--------
% image                3-D HSI
% par.blocksize        the size of each block
% par.stepsize         stepsize of each sub-cub
% -------output--------
% outpatch           switch each block into patch

% this function is the inverse of patchToimage.m
% by Wei He
%%

[M N p] = size(image);
blocksize = par.blocksize;
stepsize  = par.stepsize;

rr        = par.rr;
cc        = par.cc;
row       = par.row;
column    = par.column;
Idxx      = par.Idxx;
Numofpatch= row * column; 

outpatch = zeros(blocksize.^2,p,Numofpatch);

for idx = 1:Numofpatch
    [rowidx,columnidx] = find(Idxx==idx); % find the location in index image
    i = rr(rowidx); j = cc(columnidx);    % find the location in original image
    outpatch(:,:,idx) = reshape(image(i:i+blocksize-1,j:j+blocksize-1,:),blocksize.^2,p);
end   



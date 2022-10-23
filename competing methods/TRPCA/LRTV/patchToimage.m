function [ image] = patchToimage(outpatch,par )
% this function is to reshape the patches into 3-D image 
% -------output--------
% outpatch             3-D file, store patches 
% par.blocksize        the size of each block
% par.stepsize         stepsize of each sub-cub
% -------output--------
% image                3-D HSI

% this function is the inverse of imageTopatch.m
% by Wei He 
%%

% [M1 N1 p1] = size(outpatch);
blocksize  = par.blocksize;
stepsize   = par.stepsize;
imagesize  = par.imagesize;
p          = imagesize(3);

rr        = par.rr;
cc        = par.cc;
row       = par.row;
column    = par.column;
Idxx      = par.Idxx;
Numofpatch= row * column; 

image    = zeros(imagesize);
weight   = zeros(imagesize);

for idx = 1:Numofpatch
    [rowidx,columnidx] = find(Idxx==idx); % find the location in index image
    i = rr(rowidx); j = cc(columnidx);    % find the location in original image
    image(i:i+blocksize-1,j:j+blocksize-1,:) = image(i:i+blocksize-1,j:j+blocksize-1,:)+reshape(outpatch(:,:,idx),blocksize,blocksize,p);
    weight(i:i+blocksize-1,j:j+blocksize-1,:)= weight(i:i+blocksize-1,j:j+blocksize-1,:)+ 1;
end  

  if min(weight(:))==0
        error('the stepsize should be smaller than blocksize!');
  else
   image = image./weight;
  end
end




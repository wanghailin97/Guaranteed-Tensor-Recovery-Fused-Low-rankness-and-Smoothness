function im = sta(im)
im = double(im);
if length(size(im)) == 2
    ma = max(im(:));
    mi = min(im(:));
    im = (im - mi)/(ma - mi);
else
    ma = max(max(im,[],1),[],2);
    mi = min(min(im,[],1),[],2);
    im = (im - mi)./(ma - mi);
end
end
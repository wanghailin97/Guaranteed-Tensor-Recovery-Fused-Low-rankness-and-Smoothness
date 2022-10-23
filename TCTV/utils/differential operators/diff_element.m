function Eny = diff_element(dim,direction)
d = length(dim);
e = ones(1,d);
element1 = ones(e);
element2 = -1*ones(e);
element  = cat(direction, element1, element2);
Eny = ( abs(psf2otf(element, dim)) ).^2  ;
end

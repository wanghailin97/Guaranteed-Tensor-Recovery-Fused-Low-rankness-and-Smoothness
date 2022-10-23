function Eny = diff_element_order2(dim,direction)
d = length(dim);
e = ones(1,d);
element1 = ones(e);
element2 = -2*ones(e);
element3 = 1*ones(e);
element  = cat(direction, element1, element2, element3);
Eny = ( abs(psf2otf(element, dim)) ).^2  ;
end

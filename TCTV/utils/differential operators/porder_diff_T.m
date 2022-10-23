function DX_T = porder_diff_T(X, direction)
%compute the difference tensor (gradient map) along cerrtain direction
dim = size(X);
index_first = repmat({':'},1,ndims(X));
index_first(direction) = {1};
index_end = repmat({':'},1,ndims(X));
index_end(direction) = {dim(direction)};

slice = X(index_first{:}) - X(index_end{:});
DX  = diff(X,1,direction);
DX_T  = cat(direction,-slice,-DX);


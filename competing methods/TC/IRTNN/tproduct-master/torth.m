function orthX = torth(X)

% orth(X): orthX'*orthX=I
[n1,n2,n3] = size(X);

X = fft(X,[],3);
orthX = zeros(n1,n2,n3);
% first frontal slice
orthX(:,:,1) = gs(X(:,:,1));
% i=2,...,halfn3
halfn3 = round(n3/2);
for i = 2 : halfn3
    orthX(:,:,i) = gs(X(:,:,i));
    orthX(:,:,n3+2-i) = conj(orthX(:,:,i));
end
% if n3 is even
if mod(n3,2) == 0
    i = halfn3+1;
    orthX(:,:,i) = gs(X(:,:,i));
end
orthX = ifft(orthX,[],3);

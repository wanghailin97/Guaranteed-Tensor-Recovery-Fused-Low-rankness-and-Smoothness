%-----Low-tubal-rank Tensor Recovery from Compressed Sampling via 2D Correlted Total Viariation-----%

function [X,obj,err,iter] = LRTR_CS_2DCTV(A,b,Xsize,opts)

%% parameter setting
tol       = 1e-8; 
max_iter  = 500;
rho       = 1.1;
mu        = 1e-4;
max_mu    = 1e10;
detail    = 0;

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol      = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;         end
if isfield(opts, 'rho');         rho      = opts.rho;              end
if isfield(opts, 'mu');          mu       = opts.mu;               end
if isfield(opts, 'max_mu');      max_mu   = opts.max_mu;           end
if isfield(opts, 'detail');      detail   = opts.detail;           end

%% variables initialization
n = prod(Xsize);
I = eye(n);
invA = (A'*A+I)\I;

X   = randn(Xsize);
Z   = X;
G1  = diff_1(X);
G2  = diff_2(X);
M1  = zeros(Xsize);
M2  = zeros(Xsize);
M3  = zeros(Xsize);
M4  = zeros(n, 1);

%% FFT setting
h = Xsize(1);
w = Xsize(2);
d = Xsize(3);

Eny_x = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ; 
Eny_y = ( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;

T = Eny_x + Eny_y; 

%% main loop
iter = 0;
while iter<max_iter
    iter = iter + 1;  
    Xk = X;
    Zk = Z;
    %% Update X -- solve TV by FFT 
    H = diff_1T(mu*G1-M1)+diff_2T(mu*G2-M2);
    X = real( ifftn( fftn( mu*Z+M3+H )./(mu*(1+T)) ) );
  
    %% Updata G1 G2-- proximal operator of NN
    [G1M,nn_G1] = prox_nuclear(reshape(diff_1(X)+M1/mu,[h*w,d]),1/mu);
     G1         = reshape(G1M,[h,w,d]);
    [G2M,nn_G2] = prox_nuclear(reshape(diff_2(X)+M2/mu,[h*w,d]),1/mu);
     G2         = reshape(G2M,[h,w,d]);
    
    %% Update Z 
    vecZ = mu*invA*(mu*(X(:)+A'*b)+M3(:)+A'*M4);
    Z = reshape(vecZ, Xsize);
    
    %% Stop criterion
    dy1   = X-Z;  
    dy2   = y-A*X(:);
    chgX  = max(abs(Xk(:)-X(:)));
    dif   = norm(X(:)-Xk(:))/norm(X(:));
    chgZ  = max(abs(Zk(:)-Z(:)));
    chg   = max([chgX chgZ max(abs(dy1)) max(abs(dy2)) dif]);
    if chg < tol
        break;
    end 
    
    %% Update detail display
    if detail
        if iter == 1 || mod(iter, 10) == 0
            obj = nn_G1+nn_G2;
            err = norm(dy2,'fro');
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)]); 
        end
    end
    
    %% Update M1 M2 M3 M4 mu
    M1 = M1+mu*(diff_1(X)-G1);
    M2 = M2+mu*(diff_2(X)-G2);
    M3 = M3+mu*dy1;
    M4 = M4+mu*dy2;
    mu = min(rho*mu,max_mu);    
    
end
end
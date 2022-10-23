%-----Low-tubal-rank Tensor Completion via 2D Correlted Total Viariation-----%
function [X,obj,err,iter] =LRTC_TNN_2DCTV(M,Omega,opts) 

%% paremeters setting
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
dim = size(M);
X   = randn(dim);
G1  = diff_1(X);
G2  = diff_2(X);
E   = zeros(dim);
M1  = zeros(dim);
M2  = zeros(dim);
M3  = zeros(dim);

%% FFT setting
h = dim(1);
w = dim(2);
d = dim(3);

Eny_x = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ; 
Eny_y = ( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;

T = Eny_x + Eny_y; 

%% main loop
iter = 0;
while iter<max_iter
    iter = iter + 1;  
    Xk = X;
    Ek = E;
    %% Update X -- solve TV by FFT
    H = diff_1T(mu*G1-M1)+diff_2T(mu*G2-M2);
    X = real( ifftn( fftn( mu*(M-E)+M3+H )./(mu*(1+T)) ) );
  
    %% Updata G1 G2-- proximal operator of TNN
    [G1,tnn_G1] = prox_tnn(diff_1(X)+M1/mu,1/mu);
    [G2,tnn_G2] = prox_tnn(diff_2(X)+M2/mu,1/mu);
    
    %% Update E 
    E          = M-X+M3/mu;
    E(Omega)   = 0;
    
    %% Stop criterion
    dY   = M-X-E;    
    chgX = max(abs(Xk(:)-X(:)));
    chgE = max(abs(Ek(:)-E(:)));
    chg  = max([chgX chgE max(abs(dY(:)))]);
    if chg < tol
        break;
    end 
    
    %% Update detail display
    if detail
        if iter == 1 || mod(iter, 10) == 0
            obj = tnn_G1+tnn_G2;
            err = norm(dY(:),'fro');
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)]); 
        end
    end
    
    %% Update M1 M2 M3 mu
    M1 = M1+mu*(diff_1(X)-G1);
    M2 = M2+mu*(diff_2(X)-G2);
    M3 = M3+mu*(M-X-E);
    mu = min(rho*mu,max_mu);    
    
end
end


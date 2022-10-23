%-----Low-tubal-rank Tensor Completion via TNN plus 3D Total Viariation-----%
function [X,obj,err,iter] =LRTC_TNN_3DTV(M,Omega,lambda,opts) 

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
G3  = diff_3(X);
Y   = zeros(dim);
E   = zeros(dim);
M1  = zeros(dim);
M2  = zeros(dim);
M3  = zeros(dim);
M4  = zeros(dim);
M5  = zeros(dim);

%% FFT setting
h = dim(1);
w = dim(2);
d = dim(3);

Eny_x = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ; 
Eny_y = ( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;
Eny_z = ( abs(psf2otf([+1, -1], [w,d,h])) ).^2  ;
Eny_z = permute(Eny_z, [3, 1 2]);

T = Eny_x + Eny_y + Eny_z; 

%% main loop
iter = 0;
while iter<max_iter
    iter = iter + 1;  
    Xk = X;
    Ek = E;
    Yk = Y;
    %% Update X -- solve TV by FFT 
    H = diff_1T(mu*G1-M1)+diff_2T(mu*G2-M2)+diff_3T(mu*G3-M3);
    X = real( ifftn( fftn( mu*(M+Y-E)+M4+M5+H )./(mu*(2+T)) ) );
  
    %% Update Y -- proximal operator of TNN
    [Y,tnn_Y] = prox_tnn(X-M4/mu,lambda/mu);
    
    %% Updata G1 G2 -- soft-thresholding operator 
    G1 = prox_l1(diff_1(X)+M1/mu,1/mu);
    G2 = prox_l1(diff_2(X)+M2/mu,1/mu);
    G3 = prox_l1(diff_3(X)+M3/mu,1/mu);
    
    %% Update E 
    E          = M-X+M5/mu;
    E(Omega)   = 0;
    
    %% Stop criterion
    dY   = M-X-E;    
    chgX = max(abs(Xk(:)-X(:)));
    chgE = max(abs(Ek(:)-E(:)));
    chgY = max(abs(Yk(:)-Y(:)));
    chg  = max([chgX chgE chgY max(abs(dY(:)))]);
    if chg < tol
        break;
    end 
    
    %% Update detail display
    if detail
        if iter == 1 || mod(iter, 10) == 0
            obj = tnn_Y+norm(G1(:),1)+norm(G2(:),1)+norm(G3(:),1);
            err = norm(dY(:),'fro');
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)]); 
        end
    end
    
    %% Update M1 M2 M3 mu
    M1 = M1+mu*(diff_1(X)-G1);
    M2 = M2+mu*(diff_2(X)-G2);
    M3 = M3+mu*(diff_3(X)-G3);
    M4 = M4+mu*(Y-X);
    M5 = M5+mu*(M-X-E);
    mu = min(rho*mu,max_mu);    
    
end
end


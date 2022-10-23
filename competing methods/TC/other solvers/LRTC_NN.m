%-----Low-tubal-rank Tensor Completion via NN-----%
function [X,obj,err,iter] = LRTC_NN(M,omega,opts)

%% paremeters setting
tol      = 1e-8; 
max_iter = 500;
rho      = 1.1;
mu       = 1e-4;
max_mu   = 1e10;
detail   = 0;

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'detail');      detail = opts.detail;        end

%% variables initialization
dim = size(M);
h = dim(1);
w = dim(2);
d = dim(3);
X = zeros(dim);
X(omega) = M(omega);
E = zeros(dim);
Y = E;

%% main loop
for iter = 1 : max_iter
    Xk = X;
    Ek = E;
    % update X
    [XM,nn_XM] = prox_nuclear(reshape(M-E+Y/mu,[h*w,d]),1/mu);
     X         = reshape(XM,[h,w,d]);
    %[X,tnnX] = prox_tnn(-E+M+Y/mu,1/mu); 
    % update E
    E = M-X+Y/mu;
    E(omega) = 0;
 
    dY = M-X-E;    
    chgX = max(abs(Xk(:)-X(:)));
    chgE = max(abs(Ek(:)-E(:)));
    chg = max([chgX chgE max(abs(dY(:)))]);
    if detail
        if iter == 1 || mod(iter, 10) == 0
            obj = nn_XM;
            err = norm(dY(:));
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)]); 
        end
    end
    
    if chg < tol
        break;
    end 
    Y = Y + mu*dY;
    mu = min(rho*mu,max_mu);    
end
obj = nn_XM;
err = norm(dY(:));

 
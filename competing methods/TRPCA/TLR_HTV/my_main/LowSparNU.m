function [L1 S1 k eta beta] = LowSparNU(UU,Y,Omega,opts)

%%
if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'beta');        beta = opts.beta;            end
if isfield(opts, 'itemax');      itemax = opts.itemax;        end
if isfield(opts, 'gamma');       gamma = opts.gamma;          end
if isfield(opts, 'rho');         rho   = opts.rho;            end
if isfield(opts, 'tau');         tau   = opts.tau;            end
if isfield(opts, 'lambda');      lambda = opts.lambda;        end
if isfield(opts, 'theta');       theta  = opts.theta;         end
if isfield(opts, 'xi');          xi     = opts.xi;            end

sizeD           = size(Y);
Obser = Y.*Omega;
BarOmega = ones(sizeD) - Omega;

%% 

h               = sizeD(1);
w               = sizeD(2);
d               = sizeD(3);
%% 
Eny_x   = ( abs(psf2otf([+1; -1], [h,w,d])) ).^2  ;
Eny_y   = ( abs(psf2otf([+1, -1], [h,w,d])) ).^2  ;
Eny_z   = ( abs(psf2otf([+1, -1], [w,d,h])) ).^2  ;
Eny_z   =  permute(Eny_z, [3, 1 2]);
determ  =  Eny_x + Eny_y + Eny_z;
FD3 = beta*determ + beta;
%%
n = prod(sizeD);

%% Initial 
L0  = zeros(sizeD);
S0  = L0;
M0  = L0;
T10 = L0;
T20 = zeros(3*n,1);     % multiplier of F-D(M)
T30 = L0;



for k = 1:itemax
    %% updata Z^{k+1/2}
    Z12 = (beta*(Obser - L0 -S0) - T10).*Omega/(beta + rho) + (-L0 - S0 - T10/beta).*BarOmega;
    
    %% update L^{k+1}, F^{k+1}
    L1 = prox_utnn(UU,(Obser + M0 - S0 - Z12 -(T10 + T30)/beta)/2,1/(2*beta));
    
    diff_M     = diff3(M0(:), sizeD); 
    F1 = prox_l1(diff_M - T20(:)/beta,gamma/beta);
    
    %% update Z^{k+1}
    Z1 = (beta*(Obser - L1 -S0) - T10).*Omega/(beta + rho) + (-L1 - S0 - T10/beta).*BarOmega;
    
    %% update S^{k+1}, M^{k+1}
    S1 = prox_l1(Obser - L1 - Z1 - T10/beta,lambda/beta);
    
    diffT_p  = diffT3( T20(:) + beta*F1, sizeD );
    numer1   = reshape( diffT_p + beta*L1(:) + T30(:), sizeD);
    z        = real( ifftn( fftn(numer1) ./ (FD3) ) );
    M1       = reshape(z,sizeD);
    
   %% update T1 T2 T3
   diff_M1     = diff3(M1(:), sizeD);
   
   T1 = T10 + tau * beta * (L1 + S1 + Z1 - Obser);
   T2 = T20 + tau * beta * (F1 - diff_M1);
   T3 = T30 + tau * beta * (L1 - M1);
   
   %% stop criterion
   if k >= itemax
   etal = L1 - prox_utnn(UU,L1 - T1 - T3,1);
   etal = norm(etal(:))/(1 + norm(L1(:)) + norm(T1(:)) + norm(T3(:)));
   
   etaf = F1 - prox_l1(F1 - T2,gamma);
   etaf = norm(etaf(:))/(1+ norm(F1(:)) + norm(T2(:)));
   
   etas = S1 - prox_l1(S1 - T1,lambda);
   etas = norm(etas(:))/(1 + norm(S1(:)) + norm(T1(:)));
   
   PZ = rho*Z1.*Omega;
   etaz = PZ + T1;
   etaz = norm(etaz(:))/(1 + norm(PZ(:)) + norm(T1(:)));
   
   diffT_pst  = diffT3( T2(:), sizeD );
   numer1st   = reshape( diffT_pst + T3(:), sizeD);
   etam = norm(numer1st(:))/(1 + norm(diffT_pst(:)) + norm(T3(:)));
   
   etat1 = L1 + S1 + Z1 - Obser;
   etat1 = norm(etat1(:))/(1 + norm(Obser(:)));
   
   etat2 = F1 - diff_M1;
   etat2 = norm(etat2(:))/(1 + norm(F1(:)) + norm(diff_M1(:)));
   
   etat3 = norm(L1(:) - M1(:))/(1 + norm(L1(:)) + norm(M1(:)));
   
   etap = max([etal, etaf, etas, etaz, etam]); 
   etad = max([etat1, etat2, etat3]);
   eta = [etap, etad];
   
   if max(eta) <= tol
       break;
   end
   
   
   %% update beta
   
   kappa = etap/etad;
   if kappa > theta
       beta = xi*beta;
   elseif 1/kappa > theta
       beta = beta/xi;
   end
   
   end
   
   
   %% update variable
   T10 = T1;
   T20 = T2;
   T30 = T3;
   L0 = L1;
   S0 = S1;
   M0 = M1;
   
end

end
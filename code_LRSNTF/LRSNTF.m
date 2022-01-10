function [X, A, B, S, Out] = OurMethod(Y, opts)

max_it   = opts.max_it;  
Bmax_it  = opts.Bmax_it; 
tol      = opts.tol;              
R        = opts.R;
rho      = opts.rho;
tau      = opts.tau;
lambda   = opts.lambda;
beta     = opts.beta;
mu       = opts.mu;


Out.Res=[]; Out.PSNR=[];

%% Initiation
Nway     = size(Y);%eg 150*150*31
NwayB    = [Nway(1),Nway(2),R];
A        = abs(rand(Nway(3), R));
B        = abs(rand(NwayB));
Y3       = Unfold(Y,Nway,3);
X        = Y;
S        = zeros(Nway);
N        = zeros(Nway);
Phi       = 0;
%% Difference operator 
D= cell(1,3);     
for i = 1:3
diaga = ones(Nway(i),1);  diagb = ones(Nway(i)-1,1);
D{i}  = diag(-diaga)+diag(diagb,1);
D{i}(end,1) = 1;
end

%% D1 and D2
d1         = zeros(Nway(1),Nway(2));%eg 150*150
d1(end,1)  = 1;  d1(1,1) = -1;
d2         = zeros(Nway(1),Nway(2));
d2(1,end)  = 1;  d2(1,1) = -1;%第一行。。。
d1eig      = fft2(d1);
d2eig      = fft2(d2);
SigD12tD12 = beta*((abs(d1eig)).^2+(abs(d2eig)).^2);
SigD12tD12 = SigD12tD12(:);

for k=1: max_it
     Ak = A;   Bk = B;  Xk = X;  Sk=S; Nk=N;
     
    %% update A
    Bk3       = Unfold(Bk,NwayB,3);
    Sk3       = Unfold(Sk,Nway,3);  
    Nk3       = Unfold(Nk,Nway,3);  
    
    %%

    G         = (Y3-Sk3-Nk3)*Bk3'+ rho*Ak;
    Nn = Bk3*Bk3' + (rho + mu)*eye(size(Bk3,1));
    A = G*inv(Nn);

    
    
    %% update B
    B = GroupTV_B(Y, A, Bk, Sk, Nk, D, NwayB, SigD12tD12, rho, tau, beta, Bmax_it);
    X = double(ttm(tensor(B),A,3));%ttm就是张量乘法

    
    
    %% update S
    W = 1/(abs((Y-X-Nk+rho*Sk)/(1+rho))+eps);
    S = prox_l1((Y-X-Nk+rho*Sk)/(1+rho),W*mu/(1+rho));
    %% Update N
    N = (Y - X - S) / (1 + lambda);
   
    %% 
    Res = norm(X(:)-Xk(:))/norm(Xk(:));
    Out.Res = [Out.Res,Res];
    
    if isfield(opts, 'Xtrue')
        XT=opts.Xtrue;
        psnr = PSNR3D(XT * 255, X * 255);
        Out.PSNR = [Out.PSNR,psnr];
    end
    
    if mod(k, 2) == 0
         if isfield(opts, 'Xtrue')
            fprintf('LRSNTF: iter = %d   PSNR = %f   res = %f \n', k, psnr, Res);
         else
            fprintf('LRSNTF: iter = %d   res = %f\n', k, Res);
         end  
    end
    
    %% check stopping criterion
    if Res<tol
        break
    end
    
    
end

end



function B = GroupTV_B(Y, A, Bk, Sk, Nk, D, NwayB, SigD12tD12, rho, tau, beta, Bmax_it)

% auxiliary variable
Z1 = zeros(NwayB);
Z2 = zeros(NwayB);
% multiplier
P1=zeros(NwayB);
P2=zeros(NwayB);

for i=1:Bmax_it
    %% B subproblem
    [U2,S2,~]  = svd(A'*A);
    Sig2       = diag(S2);
    Sigg       = repmat(Sig2',NwayB(1)*NwayB(2),1)+repmat(SigD12tD12,1,NwayB(3))+rho;
    Sigg       = 1./Sigg;
    K          = double(ttm(tensor(Y-Sk-Nk),A',3))+beta*(double(ttm(tensor(Z1-P1/beta),D{1}',1))+double(ttm(tensor(Z2-P2/beta),D{2}',2)))+rho*Bk;
    K3         = Unfold(K,NwayB,3);
    temp       = Sigg.*(calF(K3',NwayB(1),NwayB(2))*U2);%calF应该是张量的DFT
    B3t        = real(calFt(temp,NwayB(1),NwayB(2)))*U2';
    B          = Fold(B3t',NwayB,3);%??起来
    
    %% Z subproblem
    Z1         = Thres_21(double(ttm(tensor(B),D{1},1))+P1/beta, tau/beta);
    Z2         = Thres_21(double(ttm(tensor(B),D{2},2))+P2/beta, tau/beta);
    %% updating P
    P1 =P1+beta*(double(ttm(tensor(B),D{1},1))-Z1);
    P2 =P2+beta*(double(ttm(tensor(B),D{2},2))-Z2);
end
end



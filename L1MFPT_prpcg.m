function bestx = L1MFPT_prpcg(M,sigma,h,c,N)
tol = 1e-3;
M = M + 1e-8;
n = size(M,1);
x0 = zeros(n,1);
g0 = grad(M,sigma,h,c,N,x0);
ite = 0;
d0 = -g0;
dMd = d0'*M*d0;
if dMd == 0
    bestx = x0;
    return;
end
b = (d0'*M*x0+sigma*h'*d0)/dMd;
k = 2*c/dMd; % c
q = N*x0;
f = N*d0;
a0 = L1MFPT_Smart(b,k,q(f~=0,:),f(f~=0,:));
g0_1 = g0;
x1 = x0+a0*d0;
if norm(x1-x0) < tol
    bestx = x1;
    return;
end
g0 = grad(M,sigma,h,c,N,x1);
x0 = x1;
while ite < 500
    beta = g0'*(g0-g0_1)/(g0_1'*g0_1);
    d0 = -g0+beta*d0;
    dMd = d0'*M*d0;
    if dMd < 1e-6
        bestx = x0;
        return;
    end
    b = (d0'*M*x0+sigma*h'*d0)/dMd;
    k = 2*c/dMd;
    q = N*x0;
    f = N*d0;
    a0 = L1MFPT_Smart(b,k,q(f~=0,:),f(f~=0,:));
    g0_1 = g0;
    x1 = x0+a0*d0;
    if norm(x1-x0) < tol
        break;
    end
    g0 = grad(M,sigma,h,c,N,x1);
    x0 = x1;
    ite = ite+1;
end
bestx = x0;
end

function right = grad(M,sigma,h,c,N,x)
% grad(x) = Mx+sigma*h-cN'sign((e-Nx)_+),where sign(0) = 0.
tol = 1e-6;
right = ones(size(N,1),1);
Nx = N*x;
Nxe = 1-Nx;
right(Nxe<tol,1) = 0;
right = M*x+sigma*h-c*N'*right;
end
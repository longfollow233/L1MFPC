function wi = L1MFPT_w(Ai,Bi,c,sigma)
X = [Ai,Bi];
XX = X * X';
mean_A = mean(Ai,2);
BarAi = Ai-mean_A;
BarBi = Bi-mean_A;
BarAiT = BarAi';
BarBiT = BarBi';
AA = BarAi*BarAiT;
[V,D] = eig(AA);
[~,p] = min(diag(D));
w0 = real(V(:,p));
iter = 0; tol = 1e-3; eps = 1e-8;
while iter < 30
    iter = iter+1;
    Aw = abs(BarAiT*w0);
    wi = w0;
    if sum(Aw) <= eps    
        break;
    end
    Gi = diag(1./(Aw+eps));
    wXw = w0'*XX*w0;
    hi = sign(wXw-1)*XX*w0;
    Fi = diag(sign(BarBiT*w0));
    M = BarAi*Gi*BarAiT;
    N = Fi*BarBiT;
    wi = L1MFPT_prpcg(M,sigma,hi,c,N);
    som = norm(wi-w0);
    if som < tol
        break;
    end
    w0 = wi;
end
end
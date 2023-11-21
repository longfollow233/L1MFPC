function bestx=L1MFPT_Exaush(b,c,q,f)
% One dimension searching: min x^2+2bx+ce'(e-q-fx)_+
% Where b,c are numbers,e,q，f are vectors and e is ones.Note: f must have no 0. 
val = (1-q)./f; 
n = length(val);
k1 = zeros(n+1,1);
k2 = zeros(n+1,1);
coefl = zeros(n,1);
coefr = zeros(n,1);
hl = zeros(n,1);
hr = zeros(n,1);
for i = 1:n
    if f(i) > 0
        coefl(i) = -f(i);
        coefr(i) = 0;
        hl(i) = 1-q(i);
        hr(i) = 0;
    elseif f(i) < 0
        coefl(i) = 0;
        coefr(i) = -f(i);
        hl(i) = 0;
        hr(i) = 1-q(i);
    end
end
[val,ind] = sort(val);
coefl = coefl(ind,1);
coefr = coefr(ind,1);
hl = hl(ind,1);
hr = hr(ind,1);
for i = 1:length(val)+1
    if i == 1
        left = -inf;
        right = val(1);
    elseif i == length(val)+1
        left = val(length(val));
        right = inf;
    else
        left = val(i-1);
        right = val(i);
    end
    if i == 1
        k1(1) = sum(coefl);
        k2(1) = sum(hl);
    elseif i == n+1
        k1(n+1) = sum(coefr);
        k2(n+1) = sum(hr);
    else 
        k1(i) = k1(i-1)-coefl(i-1)+coefr(i-1);
        k2(i) = k2(i-1)-hl(i-1)+hr(i-1);
    end
    coef1 = k1(i)*c+2*b;
    coef2 = k2(i)*c;
    x = -coef1/2;
    if x < left
        x = left;
    elseif x > right
        x = right;
    end
    opt = x^2+coef1*x+coef2;
    if i == 1
        bestx = x;
        bestf = opt;
    elseif opt < bestf
        bestf = opt;
        bestx = x;
    end
end
end
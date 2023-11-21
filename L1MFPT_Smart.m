function bestx = L1MFPT_Smart(b,c,q,f)
% One dimension searching: min x^2+2bx+ce'(e-q-fx)_+
% Where b,c are numbers,e,f,q are vectors and e is ones.Note: f must have no 0. 
val = (1-q)./f;
k = 0;
n = length(f);
tol = 1e-3;
for i = 1:n
    if -b<val(i) && f(i)>0||-b>val(i) && f(i)<0
        k = k-f(i);
    end
end
if abs(k) < tol
    bestx = -b;
    return;
else
    coefl = zeros(n,1);
    coefr = zeros(n,1);
    for i = 1:n
        if f(i) > 0
            coefl(i) = -f(i);
            coefr(i) = 0;
        elseif f(i) < 0
            coefl(i) = 0;
            coefr(i) = -f(i);
        end
    end
end
region_coef = zeros(n,1);
if k < 0
    [val,ind] = sort(val);
    coefl = coefl(ind,1);
    coefr = coefr(ind,1);
    region_coef(1) = sum(coefl);
    last_point = [];
    for i = 2:length(region_coef)
        new_coef = region_coef(i-1,1)-coefl(i-1)+coefr(i-1);
        if new_coef < 0
            region_coef(i,1) = new_coef;
        else
            last_point = i-1;
            break;
        end
    end
    if isempty(last_point)
        last_point = n;
    end
    [~,right] = min(abs(region_coef(1:last_point,1)-k));
    s = -(b+c/2*region_coef(right,1));
    while 1
        if s <= val(right,1)
            bestx = s;
            return;
        elseif right == last_point
            bestx = val(right,1);
            return;
        else
            right = right+1;
        end
        s = -(b+c/2*region_coef(right,1));
        if s <= val(right-1,1)
            bestx = val(right-1,1);
            return;
        end
    end           
else
    [val,ind] = sort(val,'descend');
    coefl = coefl(ind,1);
    coefr = coefr(ind,1);
    region_coef(1) = sum(coefr);
    last_point = [];
    for i = 2:length(region_coef)
        new_coef = region_coef(i-1,1)-coefr(i-1)+coefl(i-1);
        if new_coef > 0
            region_coef(i,1) = new_coef;
        else
            last_point = i-1;
            break;
        end
    end
    if isempty(last_point)
        last_point = n;
    end
    [~,left] = min(abs(region_coef(1:last_point,1)-k));
    s = -(b+c/2*region_coef(left,1));
    while 1        
        if s >= val(left,1)
            bestx = s;
            return;
        elseif left == last_point
            bestx = val(left,1);
            return;            
        else
            left = left+1;
        end
        s = -(b+c/2*region_coef(left,1));
        if s >= val(left-1,1)
            bestx = val(left-1,1);
            return;
        end
    end
end
end
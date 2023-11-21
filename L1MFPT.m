function ty = L1MFPT(X,iY,c,sigma,mt)

X       = X'; 
py      = iY;
[n,m]   = size(X);
flag    = 0; 
alliter = 0;

while flag==0 && alliter<30
    tic;
    alliter = alliter + 1;
    ty = py;
    L  = unique(ty);
    if length(L) <= 1
        return;
    end
    num = length(L);
    CenterX = zeros(n,num);
    Dis = zeros(m,num);
    Wti = zeros(n,mt,num);
    for i = 1:num
        Ai = X(:,ty==L(i));
        Bi = X(:,ty~=L(i));
        CenterX(:,i) = mean(Ai,2);
        for t = 1:mt
            if t==1
                wt = zeros(n,1);
            elseif norm(wt)~=0
                wt = wt/norm(wt);
            end
            wt2 = wt*wt';
            Ai = Ai-wt2*Ai;
            Bi = Bi-wt2*Bi;
            wt = L1MFPT_w(Ai,Bi,c,sigma);
            Wti(:,t,i) = wt;
        end
    end
    toc;
    for i = 1:num
        M = Wti(:,:,i)'*(X-repmat(CenterX(:,i),1,m));
        for j = 1:m
            Dis(j,i) = norm(M(:,j),1);
        end
    end
    [~,py] = min(Dis,[],2);
    if getAC(ty,py) > 0.9999
        flag = 1;
    end
end
end

function ac = getAC(ty,py)
% ty:Actual label,py:Predictive label
if size(ty,1)~=size(py,1) || size(ty,2)~=size(py,2)
    error('The size of two vector is not the same in getAC.m');
end
m = size(ty,1);
tM = zeros(m*(m-1)/2,1);
pM = zeros(m*(m-1)/2,1);
num = 1;
column = 1;
for i = 1:m-1
    column = column+1;
    j = column;
    while j ~= m+1 
        if ty(i,1) == ty(j,1)
            tM(num,1) = 1;
        end
        if py(i,1) == py(j,1)
            pM(num,1) = 1;
        end
        num = num+1;
        j = j+1;
    end
end
ac = size(find(tM==pM),1)/(m*(m-1)/2);
end
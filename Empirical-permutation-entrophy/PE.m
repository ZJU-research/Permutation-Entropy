% x - time series, Tau - delay ¦Ó , d - order of ordinal pattern d ,
% WS - size M of a sliding window
% ePE - values of the emprical permutation entrophy
function ePE = PE(x, Tau, d, WS)
load(['table', num2str(d), '.mat']);
pTbl = eval(['table', num2str(d)]);
Length = numel(x) ;
d1 = d + 1;
dTau = d * Tau;
nPat = factorial(d1);
opd = zeros(1, nPat);
ePE = zeros(1 ,Length);
op = zeros(1, d);
prevOP = zeros(1 ,Tau);
opW = zeros(1, WS);
ancNum = nPat ./ factorial(2 : d1);
peTb1 = zeros(1, WS);
peTbl(1:WS) = -(1:WS) .* log(1:WS);
peTbl(2:WS) = (peTbl(2:WS) - peTbl(1:WS-1)) ./ WS;
for iTau = 1 : Tau
    cnt = iTau;
    op(1) = (x(dTau+iTau-Tau) >= x(dTau+iTau));
    for j = 2:d
        op(j) = sum(x((d-j) * Tau + iTau) >= x((d1-j) * Tau + iTau : Tau : dTau+iTau));
    end
    opW(cnt) = sum(op .* ancNum);
    opd(opW(cnt) + 1) = opd(opW(cnt) + 1) + 1;
    for j = dTau+Tau+iTau : Tau : WS+dTau
        cnt = cnt + Tau;
        posL = 1;
        for i = j-dTau:Tau:j-Tau
            if (x(i) >= x(j))
            posL = posL + 1;
            end
        end
    opW(cnt) = pTbl(opW(cnt-Tau)*d1 + posL);
    opd(opW(cnt) + 1) = opd(opW(cnt) + 1) + 1;  
    end
    prevOP(iTau) = opW(cnt) ;
end
ordDistNorm = opd / WS;
ePE(WS+Tau*d ) = -nansum(ordDistNorm(1:nPat) .* log(ordDistNorm(1:nPat)));

iTau = 1;
iPat = 1;
for t = WS+Tau*d+1:Length
    posL = 1;
    for j = t-dTau : Tau : t-Tau
        if (x(j) >= x(t))
            posL = posL + 1;
        end
    end
    nNew = pTbl(prevOP(iTau)*d1+posL);
    nOut = opW(iPat);
    prevOP(iTau) = nNew ;
    opW(iPat) = nNew ;
    nNew = nNew + 1;
    nOut = nOut + 1;
    if nNew ~= nOut
        opd(nNew) = opd(nNew) + 1; 
        opd(nOut) = opd(nOut) - 1;
        ePE(t) = ePE(t-1) + (peTbl(opd(nNew)) - peTbl(opd(nOut) + 1 ));
    else
        ePE(t) = ePE(t-1);
    end
    iTau = iTau + 1;
    iPat = iPat + 1;
    if (iTau > Tau) 
        iTau = 1; 
    end
    if (iPat > WS)
        iPat = 1; 
    end
end
ePE = ePE (WS+Tau*d : end );
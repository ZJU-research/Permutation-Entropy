% x ? timeseries, Tau ? delay ¦Ó , d ? order of ordinal patterns d ,
% WS ? size M of a sliding window ,
% ePE ? value s of the empirical permutation entropy
function ePE = oldPE (x, Tau, d, WS)
Length = numel(x);
d1 = d + 1;
dTau = d*Tau ;
nPat = factorial(d1);
opd = zeros(1, nPat);
ePE = zeros(1 ,Length);
op = zeros(Tau, d);
opW = zeros(1, WS);
ancNum = nPat ./ factorial(2 : d1);
for iTau = 1 : Tau 
    cnt = iTau;
    op(iTau, 1) = (x(dTau+iTau-Tau) >= x(dTau+iTau)) ;
    for k = 2 : d
        op(iTau, k) = sum(x((d-k) * Tau + iTau) >= x((d1-k) * Tau + iTau : Tau : dTau+iTau)) ;
    end
    opW(cnt) = sum(op(iTau , : ) .* ancNum) + 1;
    opd (opW(cnt)) = opd (opW(cnt)) + 1;
    for t = dTau+Tau+iTau : Tau : WS+dTau
        op(iTau, 2 : d ) = op(iTau, 1 : d-1);
        op(iTau, 1) = (x(t-Tau) >= x(t));
        for j = 2 : d
            if (x(t-j*Tau) >= x(t))
                op(iTau, j) = op(iTau ,j) + 1;
            end
        end
    opNumber = sum(op(iTau , : ) .* ancNum) + 1;
    opd(opNumber) = opd(opNumber) + 1;
    cnt = cnt + Tau;
    opW(cnt) = opNumber;
    end
end
ordDistNorm = opd /WS;
ePE(WS+Tau*d ) = -nansum (ordDistNorm(1 : nPat) .* log(ordDistNorm(1 : nPat)));
iTau = 1;
iPat = 1; 
for t = WS+dTau + 1: Length 
    op(iTau ,2 : d ) = op(iTau ,1 : d-1);
    op(iTau , 1) = (x(t-Tau ) >= x(t));
    for j = 2 : d
        if (x(t-j*Tau ) >= x(t))
            op(iTau ,j) = op(iTau, j ) + 1;
        end
    end
    nNew = sum(op(iTau , : ) .* ancNum) + 1;
    nOut = opW(iPat) ;
    opW(iPat) = nNew ;
    if nNew ~= nOut
        opd(nNew) = opd(nNew) + 1;
        opd(nOut) = opd(nOut) - 1;
        ordDistNorm = opd /WS;
        ePE(t) = -nansum(ordDistNorm(1 : nPat) .* log(ordDistNorm(1 : nPat))) ;
    else
        ePE(t) = ePE(t-1) ;
    end
    iTau = iTau + 1;
    iPat = iPat + 1;
    if (iTau > Tau) 
        iTau = 1 ; 
    end
    if (iPat > WS) 
        iPat = 1 ; 
    end
end
ePE = ePE (WS+Tau*d : end ) ;
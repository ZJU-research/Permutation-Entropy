function ePE = PEeq(x, Tau, d, WS)
load(['tableEq', num2str(d), '.mat']);
opTbl = eval(['tableEq', num2str(d)]);
L = numel(x);
dTau = d*Tau ;
nPat = 1 ;
for i = 3 : 2 : 2*d+1
    nPat = nPat * i;
end
opd = zeros(1, nPat);
ePE = zeros(1 ,L);
b = zeros(Tau ,d);
prevOP = zeros(1, Tau);
opW = zeros(1 ,WS);
ancNum = ones(1 ,d);
for j = 2 : d
    ancNum(j) = ancNum(j-1)*(2*j-1);
end
peTbl(1:WS) = -(1:WS) .* log(1:WS);
peTbl(2:WS) = (peTbl(2:WS)-peTbl(1:WS-1)) ./ WS;
for iTau = 1 : Tau
    cnt = iTau;
    mOP = zeros(1 ,d);
    t = dTau+iTau;
    for j = 1 : d
        for i = j-1:-1:0
            if (i == 0 || b(iTau, i ) == 0)
                if (x(t-j*Tau) > x(t-i*Tau))
                    mOP(j) = mOP(j) + 2;
                elseif (x(t-j*Tau) == x(t-i*Tau))
                    b(iTau, j) = 1;
                end
            end
        end
    end
    mOP(1:d) = mOP(1 : d)+b(iTau ,1 : d );
    opW(cnt) = sum(mOP .* ancNum);
    opd (opW(cnt) + 1) = opd(opW(cnt) + 1) + 1;
    cnt = cnt +Tau ;
    for t = iTau +Tau*(d + 1):Tau:WS+Tau*d
        b(iTau, 2 : d ) = b(iTau, 1 : d-1);
        b(iTau ,1) = 0 ;
        posL = 1;
        eqFlag = 0;
        for i = 1:d
            if (b(iTau, i) == 0)
                if (x(t-i*Tau) > x(t))
                posL = posL + 2;
                elseif (x(t) == x(t-i * Tau))
                    eqFlag = 1;
                    b(iTau, i) = 1;
                end
            end
        end
        posL = posL+ eqFlag;
        opW(cnt) = opTbl(opW(cnt-Tau)*(2*d + 1) + posL) ;
        opd(opW(cnt) + 1) = opd(opW(cnt) + 1) + 1;
        cnt = cnt + Tau ;
    end
    prevOP(iTau) = opW(t-dTau) ;
end
OPDnorm = opd / WS;
ePE (WS+Tau*d ) = -nansum ( OPDnorm (1 : nPat ) .* log(OPDnorm(1 : nPat)));
iTau = 1;
iOP = 1;
for t = WS+Tau*d + 1:L 
    b(iTau, 2 : d ) = b (iTau, 1 : d-1) ;
    b(iTau, 1 ) = 0;
    posL = 1;
    eqFlag = 0;
    for i = 1 : d
        if (b(iTau, i) == 0)
            if (x(t-i*Tau) > x (t))
            posL = posL + 2;
            elseif (x(t) == x(t-i*Tau))
            eqFlag = 1 ;
            b(iTau ,i) = 1;
            end
        end
    end
    posL = posL+ eqFlag; 
    nNew = opTbl(prevOP(iTau)*(2*d + 1) +posL);
    nOut = opW(iOP);
    prevOP(iTau) = nNew ;
    opW(iOP) = nNew ;
    nNew = nNew + 1;
    nOut = nOut + 1;
    if nNew ~= nOut
        opd(nNew) = opd(nNew) + 1;
        opd(nOut) = opd(nOut) - 1;
        ePE(t) = ePE(t-1) + peTbl(opd(nNew)) - peTbl(opd(nOut) + 1 ) ;
    else
        ePE(t) = ePE(t-1);
    end
    iTau = iTau + 1;
    iOP = iOP + 1;
    if (iTau > Tau ) 
        iTau = 1 ; 
    end
    if (iOP > WS) 
        iOP = 1 ; 
    end
end
ePE = ePE (WS+Tau*d : end) ;
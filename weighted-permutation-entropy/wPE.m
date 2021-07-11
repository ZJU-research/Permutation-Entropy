% x - time series, Tau - delay τ , d - order of ordinal pattern d ,
% WS - size M of a sliding window
% ePE - values of the emprical permutation entrophy
function ePE = wPE(x, Tau, d, WS)
load(['table', num2str(d), '.mat']);
pTbl = eval(['table', num2str(d)]); % the precomputed table
Length = numel(x) ; %length of the time series
d1 = d + 1;
dTau = d * Tau;
nPat = factorial(d1);% amount of ordinal patterns of order d
opd = zeros(1, nPat);% distribution of ordinal patterns
ePE = zeros(1 ,Length);%empirical permutation entropy
op = zeros(1, d);% ordinal patterns
prevOP = zeros(1 ,Tau); % previous ordinal patterns
opW = zeros(1, WS); %sliding windows
ancNum = nPat ./ factorial(2 : d1); % ancillary numbers
total_var_list = zeros(1,WS); % The variance of each ordinal patterns in the sliding window
var_x_list = zeros(Tau,d+1); % The value of a segment
cnt2 = zeros(Tau); % the remove index of each Tau
total_var_cnt = 0:Tau-1; %
display_cnt = zeros(1,Tau);
for iTau = 1 : Tau
    cnt = iTau;
    cnt3 = d1;
    op(1) = (x(dTau+iTau-Tau) >= x(dTau+iTau)); %caculate the first of ordinal patterns
    var_x_list(iTau,cnt3) = x(dTau+iTau);
    cnt3 = cnt3 - 1;
    var_x_list(iTau,cnt3) = x(dTau+iTau-Tau); % put value into various content
    for j = 2:d % caculate the rest number of the first ordinal patterns
        op(j) = sum(x((d-j) * Tau + iTau) >= x((d1-j) * Tau + iTau : Tau : dTau+iTau));
        cnt3 = cnt3 - 1 ;
        var_x_list(iTau,cnt3) = x((d-j)*Tau+iTau);
    end
    opW(cnt) = sum(op .* ancNum); % the index of the first ordinal patterns
    now_mse = var(var_x_list(iTau,:)); % the various of the first ordinal patterns
    total_var_list(mod(total_var_cnt(iTau),WS)+1) = now_mse; % add it into content
    total_var_cnt(iTau) = total_var_cnt(iTau) + Tau;
    opd(opW(cnt) + 1) = opd(opW(cnt) + 1) + now_mse; % add it into the odw
    display_cnt(iTau) = display_cnt(iTau) + 1;
    for j = dTau+Tau+iTau : Tau : WS+dTau
        cnt = cnt + Tau;
        posL = 1;
        for i = j-dTau:Tau:j-Tau
            if (x(i) >= x(j))
            posL = posL + 1;
            end
        end
        var_x_list(iTau,mod(cnt2(iTau),(d+1))+1) = x(j); % add the next element and remove the last element to caculate the sum
        now_mse = var(var_x_list(iTau,:)); %  the various of the n ordinal patterns
        total_var_list(mod(total_var_cnt(iTau),WS)+1) = now_mse;
        total_var_cnt(iTau) = total_var_cnt(iTau) + Tau;
        display_cnt(iTau) = display_cnt(iTau) + 1;
        cnt2(iTau) = cnt2(iTau) + 1;
        opW(cnt) = pTbl(opW(cnt-Tau)*d1 + posL);
        opd(opW(cnt) + 1) = opd(opW(cnt) + 1) + now_mse;  
    end
    prevOP(iTau) = opW(cnt) ;
end
ordDistNorm = opd / sum(opd);
ePE(WS+Tau*d ) = -nansum(ordDistNorm(1:nPat) .* log(ordDistNorm(1:nPat)));
min_num = 99999999999999999999999;
iTau = 1;
for i = 1:Tau % get the next caculate Tau
    if min_num > display_cnt(i)
        min_num = display_cnt(i);
        iTau = i;
    end
end
iPat = 1;
for t = WS+Tau*d+1:Length
    posL = 1;
    for j = t-dTau : Tau : t-Tau
        if (x(j) >= x(t))
            posL = posL + 1;
        end
    end
    var_x_list(iTau,mod(cnt2(iTau),(d+1))+1) = x(t); % Add new elem to list
    cnt2(iTau) = cnt2(iTau) + 1;
    nNew = pTbl(prevOP(iTau)*d1+posL); % use table get next index
    nOut = opW(iPat);
    prevOP(iTau) = nNew ; % save the index into various
    opW(iPat) = nNew ;
    nNew = nNew + 1;
    nOut = nOut + 1;
    now_var = var(var_x_list(iTau,:)); %caculate various
    old_var = total_var_list(mod(total_var_cnt(iTau),WS)+1);
    total_var_list(mod(total_var_cnt(iTau),WS)+1) = now_var;
    total_var_cnt(iTau) = total_var_cnt(iTau) + Tau;
    opd(nNew) = opd(nNew) + now_var; % update opd
    opd(nOut) = opd(nOut) - old_var;
    ordDistNorm = opd / sum(opd);
    ePE(t) = -nansum(ordDistNorm(1:nPat) .* log(ordDistNorm(1:nPat))); % caculate pe
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
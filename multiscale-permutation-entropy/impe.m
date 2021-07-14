function Out_IMPE = impe(X,m,t,Scale)
%Out_IMPE = impe(X,m,t,Scale)
% Calculate the Improved Multiscale Permutation Entropy (IMPE)
% Input:    X: time series;
%           m: order of permuation entropy
%           t: delay time of permuation entropy,
%           Scale: the scale factor
% Output:
%           Out_IMPE: multiscale permuation entropy
%           eg.��scale=3,�����Ϊ��scale�ֱ�Ϊ1.2.3��impe�����ɵ���ά����           

Out_IMPE=NaN*ones(1,Scale);

% Result for scale 1 is based on the original signal
Out_IMPE(1)=pe(X,m,t);

for i=2:Scale
    Temp_PE=NaN*ones(1,i);
    
    for ii=1:i
        Xs = Multi(X(ii:end),i);
        Temp_PE(ii) = pe(Xs,m,t);  % pe computes Permutation Entropy
    end
    
    Out_IMPE(i)=mean(Temp_PE);
end


function M_Data = Multi(Data,S)
% Generate the consecutive coarse-grained time series
% Input:    Data: time series;
%           S: the scale factor
% Output:
%           M_Data: the coarse-grained time series at the scale factor S

L = length(Data);
J = fix(L/S);
M_Data = NaN*ones(1,J);

for i=1:J
    M_Data(i) = mean(Data((i-1)*S+1:i*S));
end

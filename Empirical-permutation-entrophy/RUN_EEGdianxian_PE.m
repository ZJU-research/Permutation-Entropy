[header, recorddata] = edfread('eeg10.edf');
ePE = PE(recorddata(1, :), 2, 6, 50);
% ePE = PEeq(recorddata(1, :), 2, 6, 50);
% ePE = oldPE(recorddata(1, :), 2, 6, 50);
subplot(2,1,1)
plot(recorddata(1, :));
subplot(2,1,2)
plot(ePE);
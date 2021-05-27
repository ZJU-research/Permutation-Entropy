readEEG;
ePE = PE(EEGEyeState.T7, 2, 6, 50);
subplot(2,1,1)
plot(EEGEyeState.T7);
subplot(2,1,2)
plot(ePE);
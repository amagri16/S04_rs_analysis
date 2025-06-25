%% test example of TE spectrum
%%
% Sampling frequency
Fs = 256; 

% Total time for the time series (e.g., 10 seconds)
T = 6;

% Time vector
t = 0:1/Fs:T-1/Fs;

% Frequency of the causal component
f_causal = 5; 

% Generating the causal time series (x1) at 5 Hz
X = sin(2 * pi * f_causal * t);
X2 = sin(pi * f_causal * t);

% Generating some noise to add to the second time series
noise = 0.5 * randn(size(t));

% Generating the influenced time series (x2) which is a sum of x1 and noise
% This simulates that x2 is influenced by x1 at the frequency of 5 Hz
Y = (X + noise);

% Plot the time series
figure;

subplot(2, 1, 1);
plot(t, X);
title('Causal Time Series (5 Hz)');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

subplot(2, 1, 2);
plot(t, Y);
title('Influenced Time Series');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;

% Optionally, you can also add a legend or use different line properties to distinguish the plots


%%
TOrder = 2;
FOrder = 6;
TLag = 1;
Type = 1;

[TransferEntropyM,XTSymbolMatrix,XFSymbolMatrix,YTSymbolMatrix,YFSymbolMatrix,f] = MainFunctionforFDTES(X,Y,TOrder,FOrder,TLag,Type);
%%
figure();
% low frequency gets removed during the calculation of TE
imagesc(TransferEntropyM);
sum(sum(TransferEntropyM))
title("TOrder = "+num2str(TOrder)+" "+"FOrder = "+num2str(FOrder));

%% test phase TE
data = [X' Y' X2'];
[dPTE, PTE] = PhaseTE_MF(data, [], 'otnes');

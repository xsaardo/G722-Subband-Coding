function [x,recons] = g722(filename,mode,BER)

% input speech signal
[x,Fs] = audioread(filename);

% resample to 16KHz
x = resample(x,16000/Fs,1);
x = x(1:32000);

% Quantize to 14 bits
low = min(x);
high = max(x);
step = (high - low)/(2^14 - 1);
quantization_steps = low:step:high;
x = discretize(x,quantization_steps);

% Normalize
x = x/max(x);
x = x - mean(x);

% QMF coefficients
qmf = [0.366211E-03, -0.134277E-02, -0.134277E-02, 0.646973E-02, 0.146484E-02, -0.190430E-01, 0.390625E-02, 0.441895E-01, -0.256348E-01, -0.982666E-01 0.116089E+00, 0.473145E+00];

% Analysis filtering
h0 = [qmf fliplr(qmf)];
i = 1:length(h0);
h1 = (-1).^i.*h0;
xL = downsample(conv(x,h0),2);
xH = downsample(conv(x,h1),2);

%% plot impulse response
% figure;
% subplot(2,1,1);
% stem(h0);
% xlabel('n');
% ylabel('h_0[n]');
% title('Impulse response of H_0(z)');
% subplot(2,1,2);
% stem(h1);
% xlabel('n');
% ylabel('h_1[n]');
% title('Impulse response of H_1(z)');

%% plot original/filtered signal
% subplot(3,1,1);
% plot(x);
% ylabel('Original speech signal');
% subplot(3,1,2);
% plot(xL);
% ylabel('Low subband');
% xlim([1 8000]);
% subplot(3,1,3);
% plot(xH);
% ylabel('High subband');
% xlabel('n');
% xlim([1 8000]);

% Encode signal
lowBandCodewords = lowBandADPCMEncoder(xL,step);
highBandCodewords = highBandADPCMEncoder(xH,step);

% Add transmission noise at specified BER
corruptedSignalLow = addTransmissionNoise(lowBandCodewords,BER);
corruptedSignalHigh = addTransmissionNoise(highBandCodewords,BER);

% Correct invalid codewords
corruptedSignalLow(corruptedSignalLow <= 3) = 63;

% Decode signal
rL = lowBandADPCMDecoder(corruptedSignalLow,step,mode);
rH = highBandADPCMDecoder(corruptedSignalHigh,step);

% Synthesis filtering
f0 = h0;
f1 = -h1;
vL = conv(upsample(rL,2),f0);
vH = conv(upsample(rH,2),f1);
recons = vL + vH;
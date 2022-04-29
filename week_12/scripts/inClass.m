clear;
clc;
close all;

noiseFile = './Voice.m4a'
%backgroundNoiseFile = ''

% returns sampled data, y, 
% and a sample rate for that data, Fs.
[y, Fs] = audioread(noiseFile);


T = 1/Fs;
L = length(y);
t = (0:L-1)*T;

plot(t,y)

min(find(t > 2.22376));
min(find(t > 2.8462));
%y = y(98069: 127335);

plot(t(1:length(y)),y)

numSamples = floor(L/Fs);

% hanning window
for n=0:length(y) - 1
    w(n+1) = 0.5*(1-cos(2*pi*n./length(y)-1));
end

w = w';

% convert from time to frequency
Y  = fft(w.*y);

% create single sided frequency spectrum
P2 = abs(Y/L);
P1 = P2(1:floor(L/2) + 1);
P1(2:end - 1) = 2*P1(2:end - 1);

f = Fs*(0:L/2)/L;
loglog(f, P1)

% Averaginh
amp_ave = 0;
for ii = 1:numSamples
    amp_ave = amp_ave + y(Fs*(ii - 1) + 1:Fs*ii);
end

amp_ave = amp_ave/numSamples;
L_ave = length(amp_ave);
t_ave = (0: L_ave-1)*T;
figure
plot(t_ave, amp_ave, '-g')

clear w;


for n=0:length(amp_ave) - 1
    w(n+1) = 0.5*(1-cos(2*pi*n./length(y)-1));
end

w = w';

% create sing
Y_ave = fft(w.*amp_ave);
P2_ave = abs(Y_ave/L_ave);
P1_ave = P2_ave(1:floor(L_ave/2) + 1);
P1_ave(2: end -1) = 2*P1_ave(2:end - 1);

f_ave = Fs*(0:L_ave/2)/L_ave;

figure;
loglog(f,P1, '-b');
hold on
loglog(f_ave, P1_ave, '-r')
legend('Non-averaged', 'Averaged')
hold off
figure;
semilogx(f, 20*log10(P1/20e-6), '-b')
hold on
semilogx(f_ave, 20*log10(P1_ave/20e-6), '-r')
legend('Non-averaged', 'Averaged')

dB_ave = 20*log10(P1_ave/20e-6);
% Signal to noise ratio between noise source and background
% figure
% semilogx(f, 20*log10(P1_ave/P1_Back_ave), '-b')
% hold on
% dB_diff = 20*log10(P1_ave/P1_Back_ave);

%Octave Band frequency spectrum
[OctaveData,OctaveCenterFrequencies,Flow,Fhigh] = NarrowToNthOctave(f_ave, dB_ave,3 );
figure
semilogx(OctaveCenterFrequencies,OctaveData)

% A-weighting
R_A = 12194^2*OctaveCenterFrequencies.^4./(OctaveCenterFrequencies.^2 + ...
    20.6^2)./(OctaveCenterFrequencies.^2 + 12194^2)./sqrt((OctaveCenterFrequencies.^2 + ...
    107.7^2).*(OctaveCenterFrequencies.^2 + 737.9^2))
A = 20*log10(R_A) + 2;
OctaveData_A = OctaveData + A;
figure
semilogx(OctaveCenterFrequencies,OctaveData_A, '-r')

new = 10.^(OctaveData_A./10);
totalSPL = 10*log10(sum(new))



















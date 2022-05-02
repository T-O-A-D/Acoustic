
function [P1_ave, num_samples, f, f_ave, dB_ave] = Prob_1and2(SourceFile)

[y, Fs] = audioread(SourceFile);

%y = y(98641:127335);
T = 1/Fs;
L = length(y);
t = (0:L-1)*T;

% plot(t,y)


num_samples = floor(L/Fs);

%Hanning window
for n = 0:length(y) - 1
    w(n+1) = 0.5*(1-cos(2*pi*n/(length(y)-1)));
end
w = w';

%Convert from time to frequency
Y = fft(w.*y);

%Create single sided frequency spectrum
P2 = abs(Y/L);
P1 = P2(1:floor(L/2) + 1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:L/2)/L;

% loglog(f,P1)

%Averaging
amp_ave = 0;
for ii = 1:num_samples
    amp_ave = amp_ave + y(Fs*(ii-1) + 1:Fs*ii);
end

amp_ave = amp_ave/num_samples;
L_ave = length(amp_ave);
t_ave = (0:L_ave-1)*T;

% figure;
% plot(t_ave,amp_ave)
% title("Averaging")
% xlabel("t_ave")
% ylabel("amp_ave")
clear w;

%Hanning window to averaged sample
for n = 0:length(amp_ave) - 1
    w(n+1) = 0.5*(1-cos(2*pi*n/(length(amp_ave)-1)));
end
w = w';

%Create single sided frequency spectrum of averaged samples
Y_ave = fft(w.*amp_ave);
P2_ave = abs(Y_ave/L_ave);
P1_ave = P2_ave(1:floor(L_ave/2) + 1);
P1_ave(2:end-1) = 2*P1_ave(2:end-1);

f_ave = Fs*(0:L_ave/2)/L_ave;

% figure;
% loglog(f,P1,'b-')
% hold on;
% loglog(f_ave,P1_ave,'r-')
% legend('Non-averaged','Averaged')

figure;
semilogx(f, 20*log10(P1/20e-6),'b-')
hold on;
semilogx(f_ave, 20*log10(P1_ave/20e-6),'r-')
legend('Non-averaged','Averaged')

dB_ave = 20*log10(P1_ave/20e-6);
end




function [dB_diff] = Prob_3(f, P1_ave, P1back_ave)
%Signal to noise ratio between noise source and background
figure;
semilogx(f, 20*log10(P1_ave/P1back_ave),'b')
dB_diff = 20*log10(P1_ave/P1back_ave);
end
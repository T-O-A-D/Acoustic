
function [OctaveCenterFrequencies,OctaveData,Flow,Fhigh] = Prob_4and5(f_ave, dB_ave,f_ave_BG, dB_ave_BG)
%Octave Band frequency spectrum
[OctaveData,OctaveCenterFrequencies,Flow,Fhigh] = NarrowToNthOctave(f_ave,dB_ave,3);
figure;
semilogx(OctaveCenterFrequencies,OctaveData);
hold on
[OctaveData_BG,OctaveCenterFrequencies_BG,Flow_BG,Fhigh_BG] = NarrowToNthOctave(f_ave_BG,dB_ave_BG,3);
semilogx(OctaveCenterFrequencies_BG,OctaveData_BG);
hold off
title("1/3 Octave")
xlabel("OctaveCenterFrequencies")
ylabel("OctaveData")
legend(["noise", "Background"])
end
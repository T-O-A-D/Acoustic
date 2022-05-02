
function [OctaveData_A] = Prob_6(f_ave,dB_ave_ratio)

[OctaveData,OctaveCenterFrequencies,Flow,Fhigh] = NarrowToNthOctave(f_ave,dB_ave_ratio,3);

%A-weighting
R_A = 12194^2*OctaveCenterFrequencies.^4./(OctaveCenterFrequencies.^2+...
    20.6^2)./(OctaveCenterFrequencies.^2+12194^2)./sqrt((OctaveCenterFrequencies.^2+...
    107.7^2).*(OctaveCenterFrequencies.^2+737.9^2));
A = 20*log10(R_A) + 2;
OctaveData_A = OctaveData + A;
figure;
semilogx(OctaveCenterFrequencies,OctaveData_A);
title("A-weighted")
xlabel("OctaveCenterFrequencies")
ylabel("OctaveData_A")
end
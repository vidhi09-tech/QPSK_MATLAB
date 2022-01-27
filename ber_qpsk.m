clear all;
close all;
Es_N0_dB = [-3:20];
QPSK = erfc(sqrt(0.5*(10.^(Es_N0_dB/10))))-(1/4)*erfc(sqrt(0.5*(10.^(Es_N0_dB/10)))).^2;
semilogy(Es_N0_dB,QPSK,'b.-');
axis([-3 15 10^-5 1])
grid on;
xlabel('Es_N0_dB')
ylabel('error probability')
title('symbol error probability curve for QPSK')
gtext('VIDHI KUMARI')
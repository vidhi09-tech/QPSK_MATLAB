EBn=0:20;
Ebno=10.^(EBn/10);
pe_QPSK=0.5*erfc(sqrt(Ebno));
pe_MSK=0.5*erfc(sqrt(Ebno));
semilogy(EBn,pe_QPSK,'r*-',EBn,pe_MSK,'b>-');
gtext('VIDHI KUMARI')
legend('QPSK','MSK');
xlabel('EB/N0(DB)');
ylabel('BER');

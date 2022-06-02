function [output,Pb] = QPSK_Ricean_Simu(SNR_dB,code_data)    % finds the BER/SER for the given snr
output=zeros(1,length(code_data(:)));
numofbiterror = 0;
E = 1;                                 % energy per symbol
i = 0;
SNR = 10^(SNR_dB/10);                  % signal to noise ratio
sgma = sqrt(E/SNR)/2;                  % noise variance 

% signal mapping
s00 = [1 0];
s01 = [0 1];
s11 = [-1 0];
s10 = [0 -1];

% generation of the data source, QPSK modulation
for i=1:(length(code_data(:))/2)
dsource1 = code_data(2*i-1);dsource2 = code_data(2*i); 
    
% detection and the probability of error calculation
n = sgma*randn(size(s00));                                % 2 normal distributed r.v with 0, variance sgma
Ric1=0.316223*randn(size(s00));
Ric2=0.8945+0.316223*randn(size(s00));
Ric=sqrt(Ric1.*Ric1+Ric2.*Ric2);
if((dsource1==0)&(dsource2==0)),r = Ric.*s00+n;           % go through Ricean
elseif((dsource1==0)&(dsource2==1)),r = Ric.*s01+n;
elseif((dsource1==1)&(dsource2==0)),r = Ric.*s10+n;
else r = Ric.*s11+n;
end;

c00 = dot(r,s00);c01 = dot(r,s01);
c10 = dot(r,s10);c11 = dot(r,s11);
c_max = max([c00 c01 c10 c11]);

if (c00==c_max) decis1 = 0;decis2 = 0;
elseif(c01==c_max) decis1 = 0;decis2 = 1;
elseif(c10==c_max) decis1 = 1;decis2 = 0;
else decis1 = 1;decis2 = 1;
end;

if(decis1 ~= dsource1)
  numofbiterror = numofbiterror+1;end;
if(decis2~=dsource2)
  numofbiterror = numofbiterror+1;end;
output(2*i-1) =decis1;output(2*i) = decis2; 

end;
Pb = numofbiterror/(2*i);
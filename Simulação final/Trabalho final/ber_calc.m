

clear; clc; close all;

n=31;
k=21;
t=2;
M=16;
ebn0=[0:1:15];

% n=63;
% k=51;
% t=2;
% ebn0=20.89;
% M=8;

esn0=log2(M)*ebn0*k/n;

pe=erfc(sqrt(esn0)*sin(pi/M));

pc=pe/log(M);
BER=0;
for j=t+1:1:n;
    BER=BER+j*(factorial(n)/(factorial(j)*factorial(n-j))).*pc.^j.*(1-pc).^(n-j)/n;
end
    
plot(ebn0,2*BER);
semilogy(ebn0,BER,'LineWidth',2)
axis([0 15 1e-3 1]);
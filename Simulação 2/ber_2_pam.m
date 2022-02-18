% -------------------------------------------------------------------------
% Simulacao da BER de um sistema 2-PAM no canal AWGN
% -------------------------------------------------------------------------

clear; clc; close all;

% -------------------------------------------------------------------------
% Parametros
% -------------------------------------------------------------------------
M= 2;

Fs=44100;           % Frequencia de Amostragem do Sinal Continuo (Hz)
Ts=1/Fs;            % Periodo de Amostragem (s)

oversampling=10;    % Fator Fs/R

R=Fs/oversampling;  % Taxa de Transmissao em simbolos/s (baud rate)
T=1/R;              % Periodo de Simbolo (s)

del=25;             % Resposta do filtro formatador se estende por (2*del) per?odos de simbolo
                    % Numero de amostras do filtro formatador: N=2*(del*oversampling)+1

rolloff=0.5;        % Fator de rolloff dos filtros Tx e Rx


ebn0=[0:1:15];

for k=1:1:length(ebn0)
% -------------------------------------------------------------------------
% Geração Aleatória de Símbolos 2-PAM
%--------------------------------------------------------------------------

symbols = randi(2,60000,1);
a=[-1 1];                  %alfabeto 8-PAM
pam=a(symbols);


% -------------------------------------------------------------------------
% Filtro de Tx+Rx
% Formatacao de Pulso - Tipo Raiz Quadrada do Cosseno Elevado no Tx e Rx         
% -------------------------------------------------------------------------
filtro=rcosfir(rolloff,del,oversampling,1,'sqrt');

%Formatar e transmitir os simbolos PAM

sinal_tx=upsample(pam,oversampling);           % Realiza Upsampling
sinal_tx_filtrado=conv(sinal_tx,filtro);       % Sinal Filtrado de Transmissao

% -------------------------------------------------------------------------
% Canal Passa Baixas Simulado usando um Filtro de Butterworth
% -------------------------------------------------------------------------

fc=R/2*(1+rolloff);                 % Largura de banda do sinal transmitido
[bn,an]=butter(5,2*fc/Fs);          % Filtro passa-baixas de largura de banda fc

sinal_rx=filtfilt(bn,an,sinal_tx_filtrado);


% -------------------------------------------------------------------------
% Receptor (Filtro Casado)
% -------------------------------------------------------------------------

sinal_rx_ruido=awgn(sinal_rx,(ebn0(k)+10*log10(log2(M))-10*log10(oversampling/2)),'measured');

sinal_rx_casado=conv(sinal_rx_ruido,filtro);

pam_rx=downsample(sinal_rx_casado,oversampling);     
pam_rx=pam_rx(del*2+1:length(pam_rx)-del*2);        


%Estimacao dos simbolos PAM recebidos

alphabet=[-1 1];
symbols_rx_quant=quantalph(pam_rx,alphabet);
symbols_rx=(symbols_rx_quant+1)/2+1;

[Erros,BER] = biterr(symbols,symbols_rx)
ber(k)=BER;

EbN0=10^(ebn0(k)/10);

ber_teorica(k)=0.5*erfc(sqrt(EbN0));

end
semilogy(ebn0,ber,'b.-')
hold on
semilogy(ebn0,ber_teorica,'mx-')
axis([1 15 10^-7 1])

grid
ylabel('BER'); 
xlabel('Eb/No (dB)'); 
title('Desempenho da Modulação 2-PAM no Canal AWGN');
legend('Simulação','Teórico');
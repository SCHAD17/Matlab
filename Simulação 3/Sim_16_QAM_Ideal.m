% -------------------------------------------------------------------------
% Simulacao de um Sistema de Transmissao Digital em banda Base
% Será transmitida uma imagem usando simbolos 16-QAM
% -------------------------------------------------------------------------

clear; clc; close all;

% -------------------------------------------------------------------------
% Parametros
% -------------------------------------------------------------------------

Fs=44100;           % Frequencia de Amostragem do Sinal Continuo (Hz)
Ts=1/Fs;            % Periodo de Amostragem (s)

oversampling=10;    % Fator Fs/R

R=Fs/oversampling;  % Taxa de Transmiss?o em s?mbolos/s (baud rate)
T=1/R;              % Periodo de Simbolo (s)

del=25;             % Resposta do filtro formatador se estende por (2*del) periodos de simbolo
                    % Numero de amostras do filtro formatador: N=2*(del*oversampling)+1

rolloff=0.5;        % Fator de rolloff dos filtros Tx e Rx

snr=20;             % relação sinal ruido
% -------------------------------------------------------------------------
% Leitura da Imagem e mapeamento 4-PAM
%--------------------------------------------------------------------------

im_in=imread('shuttle_80x60.tif');  %Leitura da imagem a ser transmitida
%im_in=imread('lenna512.tif');

L=8;
[size_r,size_c]=size(im_in);
im_size=size_r*size_c;
im_vec=reshape(im_in,1,im_size);

bit_matrix=de2bi(im_vec);
bit_per_symbol=2;               %constelacao 2^bit_per_symbol-PAM, neste caso 4-PAM
bit_symbols=reshape(bit_matrix, im_size*L/bit_per_symbol, bit_per_symbol);
symbols=bi2de(bit_symbols);     %dibits a ser transmitidos
alphabet=[-3 -1 1 3];           %alfabeto 4-PAM
symbols=symbols+1;
pam=alphabet(symbols);

pam_I=pam(1:2:length(pam));   % seguencia I (en fase) e Q (em quadratura)
pam_Q=pam(2:2:length(pam));   % dos simbolos QAM

%---------------------------------------------------------------------
filtro=rcosfir(rolloff,del, oversampling, 1, 'sqrt');

%Formatar e transmitir os simbolos QAM

sinal_tx_I=upsample(pam_I, oversampling);      % Realiza Upsampling
sinal_tx_Q=upsample(pam_Q, oversampling);
sinal_tx_filtrado_I=conv(sinal_tx_I, filtro);     % Sinal Filtrado de Transmissão
sinal_tx_filtrado_Q=conv(sinal_tx_Q, filtro);

sinal_QAM=sinal_tx_filtrado_I+1j*sinal_tx_filtrado_Q;

%----------------------------------------------------------------------
% Canal Passa Baixas Simulado usando um filtro de Butterworth
%-----------------------------------------------------------------------

fc=R/2*(1+rolloff);                   %Largura de banda do sinal transmitido`
[bn, an]=butter(9, 2*fc/Fs);          % Filtro passa-baixas de largura de banda fc
sinal_rx=filtfilt(bn, an, sinal_QAM); %Sinal na saida do 'canal'

sinal_rx = awgn(sinal_rx, snr, 'measured');

figure(1);
subplot(2,2,1);

% Densidade Espectral de Potência do sinal recebido
[Pxx,F]=pwelch(sinal_QAM/max(abs(sinal_QAM)),[],[],[],Fs,'twosided');
plot((F-Fs/2)/1000,10*log10(fftshift(Pxx)));
grid; xlabel('Frequência (kHz)'); ylabel('dBm/Hz');
xlim([-(Fs/2)/1000 (Fs/2)/1000]);
title('Sinal Transmitido')

subplot(2,2,2);

% Densidade Espectral de Potência do sinal recebido
[Pxx,F]=pwelch(sinal_rx/max(abs(sinal_rx)),[],[],[],Fs,'twosided');
plot((F-Fs/2)/1000, 10*log10(fftshift(Pxx)));
grid; xlabel('Frequencia (kHz)'); ylabel('dBm/Hz');
xlim([-(Fs/4)/1000 (Fs/4)/1000]);
title('Sinal na Saida do Canal')

% -------------------------------------------------------------------------
% Receptor (Filtro Casado)
% -------------------------------------------------------------------------

sinal_rx_casado=conv(sinal_rx,filtro);       %Filtro casado

sinal_rx_casado_I=conv(real(sinal_rx),filtro); %Filtro casado em fase
sinal_rx_casado_Q=conv(imag(sinal_rx),filtro); %Filtro casado em Quadratura


subplot(2,2,3);
[H,F]=freqz(bn,an,2048,'whole',Fs);
gain=20*log10(fftshift(abs(H)));
plot((F-Fs/2)/1000,gain); grid;
axis([-Fs/(4*1000) Fs/(4*1000) -50 10]);
xlabel('Frequencia (kHz)');
ylabel('Ganho (dB)');
title('Resposta de Amplitude do Canal');


subplot(2,2,4);                    

%Densidade Espectral de Potencia do sinal recebido
%[Pxx,F]=pwelch(sinal_rx_casado/max(abs(sinal_rx)),[],[],[],Fs,'twosided');
%plot((F-Fs/2)/1000,10*log10(fftshift(Pxx))); 
%grid; xlabel('Frequencia (kHz)'); ylabel('dBm/Hz');
%xlim([-(Fs/4)/1000 (Fs/4)/1000]);
%title('Sinal na Saida do Filtro Casado');


%qam_rx=upfirdn(sinal_rx_casado,oversampling);  

qam_rx=downsample(sinal_rx_casado,oversampling); 
qam_rx=qam_rx(del*2+1:length(qam_rx)-del*2);  

pam_rx_I=downsample(sinal_rx_casado_I,oversampling);
pam_rx_Q=downsample(sinal_rx_casado_Q,oversampling);
pam_rx_I=pam_rx_I(del*2+1:length(pam_rx_I)-del*2);
pam_rx_Q=pam_rx_Q(del*2+1:length(pam_rx_Q)-del*2);


%--------------------------------
eyediagram(sinal_rx_casado(1,1000:6000),2*oversampling);
sPlotFig = scatterplot(sinal_rx_casado,10,0);
hold on


%Estimação dos simbolos PAM recebidos
symbols_rx_quant=quantalph(qam_rx,alphabet);

a=[0 1 2 3];
symbols_rx=(symbols_rx_quant+3)/2+1;
symbols_rx=a(symbols_rx);                      %Reconstrucao da imagem no receptor
bit_symbols_rx=de2bi(symbols_rx);
bit_matrix_rx=reshape(bit_symbols_rx, im_size, L);
im_vec_rx=bi2de(bit_matrix_rx);
%im_vec=reshape(im_in,1,im_size);
im_in_rx=reshape(im_vec_rx,size_r,size_c);     

figure(2);

subplot(2,1,1);        %Visualizacao da imagem transmitida
colormap(gray);
h=image(im_in);
set(h,'CDataMapping','scaled')
axis('equal');
title('Imagem Transmitida');
hold;

subplot(2,1,2);                                %Visualizacao da imagem recebida
colormap(gray);
h=image(im_in_rx);
set(h,'CDataMapping','scaled')
axis('equal');
title('Imagem Recebida');
[num_erros,BER]=symerr(bit_matrix_rx,bit_matrix)

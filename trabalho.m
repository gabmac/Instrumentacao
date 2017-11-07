close all;
clc;
clear all;
% y  = dados amostrados;
%Fs = taxa de amostragem dada em Hz

%lendo a entrada
[y,Fs] = audioread('T1E1.wav');

%função para ouvir o som

%sound(y,Fs);

%audioinfo('T1E1.wav')

%Perído de Amostragem
Ts = 1/Fs;
%tempo total
tempo_total = Ts*length(y);
%vetor tempo
t = 0:Ts:tempo_total-Ts;
%sinal
plot(t,y);
title('Sinal');
ylabel('Magnitude do Sinal')
xlabel('Tempo')
figure;

%Numero de elementos da FFT  - Escolher uma potencia de 2
N = 32768;

%FFT

Y_fft = fft(y,N);
mY_fft = abs(Y_fft);
df = 1/(N*Ts);
f = [0:df:(N-1)*df]';
%Como a transformada é unidimensional, basta plotar metade do espectro

%Magnitude
plot(f,mY_fft);
xlim([0 (N-1)*df/2])
title('Transformada de Fourier (Magnitude)');
ylabel('Magnitude |\it{Y}\rm|');
xlabel('Frequencia \it{f} \rm(Hz)')
figure;


%fase
pY_fft = angle(Y_fft);
plot(f,pY_fft);
xlim([0 (N-1)*df/2])
title('Fase da Transformada de Fourier')
xlabel('Frequencia \it{f} \rm(Hz)')
ylabel('Fase /\it{Y} \rm(rad)')
figure;

%potencia
Pyy = Y_fft.*conj(Y_fft);
plot(f,Pyy);
xlim([0 (N-1)*df/2])
title('Potência')
ylabel('Potencia Pyy')
xlabel('Frequencia \it{f} \rm(Hz)')

 
%quantidade de janelamentos no sinal
janela = 8;

%tamanho em potência de 2, que mais se aproxima do sinal quando dividimos o
%tamanho do sinal pela janela
tam_janela = 4096;
figure
hold all;

%eixo das frequências
df1 = 1/(tam_janela*Ts);
f1 = [0:df1:(tam_janela-1)*df1]';

%janelamento
for i =0:1:janela-1
    y_n = y(i*tam_janela+1:(i+1)*(tam_janela));
    y_n_fft = abs(fft(y_n,tam_janela));
    plot(f1,y_n_fft);
    xlim([0 (tam_janela-1)*df/2])
    
end   
title('Transformada de Fourier Janelada (Magnitude)');
ylabel('Magnitude |\it{Y}\rm|');
xlabel('Frequencia \it{f} \rm(Hz)')
hold off;

figure
%spectogram do sinal, não é mostrado todas as harmônicas
spectrogram(y,[],[],(0:df:1024*df),Fs,'yaxis');
colormap bone;


n = 0:Ts:0.1888;
y_nf = 0.2*sin(2*pi*174.614105*n);
n =0.1888+Ts:Ts:2*0.1888;
y_nf(length(y_nf)+1:length(n)+length(y_nf)) = 0.2*sin(2*pi*184.997208*n);
n =2*0.1888+Ts:Ts:3*0.1888;
y_nf(length(y_nf)+1:length(n)+length(y_nf)) = 0.2*sin(2*pi*195.997711*n);
n =3*0.1888+Ts:Ts:1.510;
y_nf(length(y_nf)+1:length(n)+length(y_nf)) = 0.2*sin(2*pi*207.652344*n);

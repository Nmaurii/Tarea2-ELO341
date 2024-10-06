%Tarea 2 Teoria de cominicaciones digitales.

%Nicolas Vergara Fecha 02/10/2024

clc;
clear all;
close all;

%----------------ETAPA 1 DE SAMPLEO

f = 1000; % frecuancia de senal banda base 1000 Hz
n_muestras = 1000;
fs = 2*f*n_muestras; % frecuencia de muestreo (Nyquist)

T = 1/f; % periodo de un ciclo
N = (T*fs); %cantidad de muestas

t = linspace(0,T,N);%intervalo de valores.

% senal banda base a modular sampleada
A = 1;%amplitud
y = A*sin(2*pi*f*t);



%-----------------ETAPA 2 DE CUANTIZACION

n_bits = 3;
L = 2^n_bits; %vendria siendo el L en 2mp/L que nos da la cantidad de niveles.
delta = 2*A/L;%delta es la separacin enter niveles.
e = delta/2; %epsilon es el error medio 

%se define los arrays error y vector_pcm que almacena datos de cuantizacion
error = [];
vector_pcm = [];
SNR_valores = [];
%S Y N de senal
p_senal = mean(y.^2);
p_ruido = A^2/(3*L^2);

niveles = -A + delta/2 : delta : A- delta/2; 

% Inicializar la señal cuantizada
vector_pcm = zeros(size(y));

for i = 1:length(y)

    dif = abs(y(i) - niveles); 
    [minimo, idx] = min(dif);  
    vector_pcm(i) = niveles(idx); 

    error(i) = abs(y(i) - vector_pcm(i));
    SNR_valores(i) = 10*log10((y(i)^2)/p_ruido); 

end


SNR = 10*log10(p_senal/p_ruido);
disp(y(1999));


%********Grafica

figure(1);

plot(t, y, '.-', 'DisplayName', 'Señal Sampleada'); 
hold on;
plot(t, vector_pcm, 'r--', 'DisplayName', 'Señal Cuantizada'); 
plot(t, error, '-', 'DisplayName', 'Error'); 
hold off;

xlabel('Tiempo 1[ms]');
title('Señal de 1000 Hz y Señal Cuantizada');
legend('show');
grid on;
xlim([0, T]);


%SNR de senal seno banda base.
figure(2);
plot(t, SNR*ones(1, length(t)), 'b-', 'DisplayName', 'SNR '); % Línea horizontal
hold on;
plot(t(1:1999), SNR_valores(1:1999), '-', 'DisplayName', 'SNR'); 
hold off;

xlabel('Tiempo 1[ms]');
ylabel('db');
title('SNR');
legend('show');
grid on;
xlim([0, T]);
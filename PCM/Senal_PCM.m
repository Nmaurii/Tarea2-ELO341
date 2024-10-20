% Tarea 2 Teoria de comunicaciones digitales.
% Nicolas Vergara Fecha 02/10/2024

clc;
clear all;
close all;

function [senal_cuantizada, error] = Cuantizar(senal, n_bits, mp)
    L = 2^n_bits; 
    delta = 2 * mp / L;
    e = delta / 2; 
    niveles = -mp + e : delta : mp - e; 

    for i = 1:length(senal)
        dif = abs(senal(i) - niveles); 
        [~, idx] = min(dif);  
        senal_cuantizada(i) = niveles(idx); 
        error(i) = abs(senal(i) - senal_cuantizada(i));
    end
end

%----------------ETAPA 1 DE SAMPLEO
f = 1000; % Frecuencia de la señal banda base 1000 Hz
n_muestras = 1000;
fs = 2 * f * n_muestras; % Frecuencia de muestreo (Nyquist)

T = 1 / f; % Periodo de un ciclo
N = (T * fs); % Cantidad de muestras

t = linspace(0, T, N); % Intervalo de valores
A = 1; % Amplitud
senal_moduladora = A * sin(2 * pi * f * t);

%----------------ETAPA 3 ANALISIS DE SNR
n_bits_range = 3:5; % intervalor de profundidad de bits
snr_empirica = zeros(size(n_bits_range));
snr_teorica = 6.02 * n_bits_range + 1.76; % formula para SNR teórica
mp = A; % maximo  de la señal moduladora

for idx = 1:length(n_bits_range)

    n_bits = n_bits_range(idx);
    [vector_pcm, error] = Cuantizar(senal_moduladora, n_bits, mp);
    potencia_senal = mean(senal_moduladora .^ 2);
    potencia_ruido = mean((vector_pcm - senal_moduladora) .^ 2); 
    snr_empirica(idx) = 10 * log10(potencia_senal / potencia_ruido); %snr empiroco

end
figure;
plot(n_bits_range, snr_teorica, '-o', 'LineWidth', 2, 'DisplayName', 'SNR Teórica (dB)');
hold on;
plot(n_bits_range, snr_empirica, '-x', 'LineWidth', 2, 'DisplayName', 'SNR Empírica (dB)');
xlabel('Profundidad de Bits');
ylabel('SNR (dB)');
title('SNR Teórica y Empírica para profundidades de bits de 3 a 5');
legend('show');
grid on;


%Tarea 2 Teoria de cominicaciones digitales.

%Nicolas Vergara Fecha 04/10/2024

clc;
clear all;
close all;

%utilizando la misma senal banda base del item PCM    
function d_q = DPCM(m,orden,a_n)
    k = orden+1;
    m_q = zeros(1,length(m));
    d_q = zeros(1,length(m));

    for i = k:length(m)

        prediccion = a_n(1)*m_q(i-1) + a_n(2)*m_q(i-2) + a_n(3)*m_q(i-3);
        d_q(i) = m(i) - prediccion;
        m_q(i) = d_q(i) - prediccion;
    end
end

%
function [q,error,SNR] = Cuantizar(y,n_bits,mp)
    
    L = 2^n_bits; 
    delta = 2*mp/L;
    e = delta/2; 
    niveles = -mp + e : delta : mp - e; 

    for i = 1:length(y)
        dif = abs(y(i) - niveles); 
        [minimo, idx] = min(dif);  
        q(i) = niveles(idx); 

        SNR(i) = 10*log10((y(i)^2)/(mp^2/(3*(n_bits)^2))); 
        error(i) = abs(y(i) - q(i));
    end
end


% function [q,error] = Cuantizar(y,n_bits,mp)
% 
%     L = 2^n_bits; 
%     delta = 2*mp/L;
%     e = delta/2; 
%     niveles = -mp + e : delta : mp - e; 
% 
%     for i = 1:length(y)
%         dif = abs(y(i) - niveles); 
%         [minimo, idx] = min(dif);  
%         q(i) = niveles(idx); 
% 
%         error(i) = abs(y(i) - q(i));
%     end
% end


%demodulacion DPCM 
function m = DPCD(d_q, orden, a_n)
    k = orden + 1;
    m = zeros(1, length(d_q));

    for i = k:length(d_q)
        prediccion = a_n(1) * d_q(i-1) + a_n(2) * d_q(i-2) + a_n(3) * d_q(i-3);
        m(i) = d_q(i) + prediccion; 
    end
end

%****** TRANSMISOR

%----------SAMPLEO

f = 1000; % frecuancia de senal banda base 1000 Hz
n_muestras = 1000;
fs = 2*f*n_muestras; % frecuencia de muestreo (Nyquist)

T = 1/f; % periodo de un ciclo
N = (T*fs); %cantidad de muestas

t = linspace(0,T,N);%intervalo de valores.

% senal a modular
A = 1;%amplitud
m = A*sin(2*pi*f*t);

%---------MODULACION DPCM
senalDPCM_plana = DPCM(m,3,[0.7071,0.5,0.25]);

%----------CUANTIZACION
mp = max(senalDPCM_plana); %dado que DPCM contrae la senal, se recalcula el maximo de la senal
disp(mp);
cantidad_bits = 4;

%la senalDPCM_cuantizada corresponde a la TX
[senalDPCM_cuantizada,error,SNR] = Cuantizar(senalDPCM_plana,cantidad_bits,mp);

%****** RECEPTOR

senal_demodulada = DPCD(senalDPCM_cuantizada,3,[0.7071,0.5,0.25]);

figure;

plot(t, m, '.-', 'DisplayName', 'Señal Sampleada'); 
hold on;
plot(t, senalDPCM_plana, '.', 'DisplayName', 'Señal DPCM Plana TX'); 
plot(t, senalDPCM_cuantizada, 'r--', 'DisplayName', 'Señal DPCM cuantizada TX'); 
hold off;

xlabel('Tiempo 1[ms]');
title('Señal Tx');
legend('show');
grid on;
xlim([0, T]);

figure;

plot(t, senal_demodulada, 'r--', 'DisplayName', 'Señal DPCM demodulada RX'); 
hold on;


xlabel('Tiempo 1[ms]');
title('Señal Demodulada Rx');
legend('show');
grid on;
xlim([0, T]);

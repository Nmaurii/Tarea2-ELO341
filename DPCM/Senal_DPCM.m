%Tarea 2 Teoria de cominicaciones digitales.

%Nicolas Vergara Fecha 04/10/2024

clc;
clear all;
close all;

function [valor_cuantizado, error] = Cuantizar(valor, n_bits, mp)

    L = 2^n_bits; %n umero de niveles
    delta = 2 * mp / L; % tamano de cada nivel
    e = delta / 2; % error medio de los delta
    niveles = -mp + e : delta : mp - e; 
    dif = abs(valor - niveles); 
    [~, idx] = min(dif); % encuentra el indice del nivel mas cercano
    valor_cuantizado = niveles(idx); % valor cuantizado
    error = abs(valor - valor_cuantizado); % error de cuantización
end

%Funcion modulacion DPCM
function d_q = DPCM(m, orden, a_n, n_bits, mp)
    k = orden + 1;
    m_q = zeros(1, length(m));
    d_q = zeros(1, length(m));

    for i = k:length(m)
        prediccion = a_n(1) * m_q(i-1) + a_n(2) * m_q(i-2) + a_n(3) * m_q(i-3);
        d_q(i) = Cuantizar(m(i) - prediccion, n_bits, mp);
        m_q(i) = prediccion + d_q(i); % reconstrucción correcta en el transmisor
    end
end

% Funcion demodulacion DPCM
function m = DPCD(d_q, orden, a_n)
    k = orden + 1;
    m = zeros(1, length(d_q));
    m_q = zeros(1, length(d_q)); 

    for i = k:length(d_q)
        % predicción usando los coeficientes a_n
        prediccion = a_n(1) * m_q(i-1) + a_n(2) * m_q(i-2) + a_n(3) * m_q(i-3);
        m_q(i) = prediccion + d_q(i); % reconstrucción en el receptor
        m(i) = m_q(i); % senal demodulada
    end
end

%****** TRANSMISOR

%----------SAMPLEO
f = 1000; % Frecuencia de señal banda base 1000 Hz
n_muestras = 1000;
fs = 2*f*n_muestras; % Frecuencia de muestreo para 1000 muestras por ciclo

T = 1 / f; % Periodo de un ciclo
N = round(T * fs); % Cantidad de muestras

t = linspace(0, T, N); % Intervalo de valores.

% Senal moduladora
A = 1; % Amplitud
moduladora = A * sin(2 * pi * f * t);

%----------CUANTIZACION
cantidad_bits = 3;

%---------MODULACION DPCM
senalDPCM = DPCM(moduladora, 3, [0.7071, 0.5, 0.25], cantidad_bits, A);

%*******RECEPTOR
senal_demodulada = DPCD(senalDPCM, 3, [0.7071, 0.5, 0.25]);

%********Grafica

figure;
plot(t, moduladora, 'r-', 'DisplayName','Señal m[k]'); 
hold on;
plot(t, senalDPCM, 'b-', 'DisplayName', 'Señal dq[k]'); 

plot(t, senal_demodulada, 'y-', 'DisplayName', 'Señal mqr[k]');
hold off;


xlabel('tiempo [s]');
ylabel('Amplitud [V]');
title('Señal DPCM con profundidad de 3 bits');
legend('show');
grid on;
xlim([0, T]);
ylim([-1, 1]); 

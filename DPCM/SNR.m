
function [valor_cuantizado, error] = Cuantizardor(valor, n_bits, mp)
    L = 2^n_bits; % número de niveles
    delta = 2 * mp / L; % tamaño de cada nivel
    e = delta / 2; % error medio de los delta
    niveles = -mp + e : delta : mp - e; 
    dif = abs(valor - niveles); 
    [~, idx] = min(dif); % encuentra el índice del nivel más cercano
    valor_cuantizado = niveles(idx); % valor cuantizado
    error = abs(valor - valor_cuantizado); % error de cuantización
end

% Función modulación DPCM
function d_q = DPCM(m, orden, a_n, n_bits, mp)
    k = orden + 1;
    m_q = zeros(1, length(m));
    d_q = zeros(1, length(m));

    for i = k:length(m)
        prediccion = a_n(1) * m_q(i-1) + a_n(2) * m_q(i-2) + a_n(3) * m_q(i-3);
        d_q(i) = Cuantizardor(m(i) - prediccion, n_bits, mp);
        m_q(i) = prediccion + d_q(i); % reconstrucción correcta en el transmisor
    end
end

% Función demodulación DPCM
function m = DPCMD(d_q, orden, a_n)
    k = orden + 1;
    m = zeros(1, length(d_q));
    m_q = zeros(1, length(d_q)); 

    for i = k:length(d_q)
        % predicción usando los coeficientes a_n
        prediccion = a_n(1) * m_q(i-1) + a_n(2) * m_q(i-2) + a_n(3) * m_q(i-3);
        m_q(i) = prediccion + d_q(i); % reconstrucción en el receptor
        m(i) = m_q(i); % señal demodulada
    end
end

%****** TRANSMISOR

%----------SAMPLEO
f = 1000; % Frecuencia de señal banda base 1000 Hz
n_muestras = 1000;
fs = 2*n_muestras*f; % Frecuencia de muestreo para 1000 muestras por ciclo
T = 1 / f; % Periodo de un ciclo
N = round(T * fs); % Cantidad de muestras
t = linspace(0, T, N); % Intervalo de valores.

% Señal moduladora
A = 1; % Amplitud
moduladora = A * sin(2 * pi * f * t);

% Valores de bits de 3 a 8
bit_values = 3:8;
SNR_teorico = zeros(1, length(bit_values));
SNR_empirico = zeros(1, length(n_bits));

for idx = 1:length(bit_values)

    cantidad_bits = bit_values(idx);
    %--------- DPCM
    senalDPCM = DPCM(moduladora, 3, [0.7071, 0.5, 0.25], cantidad_bits, A);

    %*******RECEPTOR
    senal_demodulada = DPCMD(senalDPCM, 3, [0.7071, 0.5, 0.25]);

    %******* SNR
    % Error de cuantizacion para DPCM
    error_cuantizacion = abs(moduladora - senal_demodulada); % Error entre la señal original y la demodulada
    sigma_Q_DPCM = var(error_cuantizacion); % Varianza del error de cuantización

    % Varianza de la senal original
    sigma_X = var(moduladora); % Varianza de la señal original

    % Varianza del error de predicción
    error_prediccion = moduladora - senal_demodulada; % Error de prediccion
    sigma_E = var(error_prediccion); % Varianza del error de predicción
    sigma_D = var(senalDPCM); % varianza de la señal DPCM

    % SNR
    ganancia= sigma_X^2 / sigma_Q_DPCM; % ganancia de predicción
    SNR_Q = (sigma_X^2 / sigma_Q_DPCM) * (sigma_Q_DPCM / sigma_D); % SNR total

    SNR_empirico(idx) = 10 * log10(SNR_Q * ganancia); % Guardar en dB
    SNR_teorico(idx) = 6.02 * cantidad_bits + 10 * log10(sigma_X / sigma_D) + 1.76; % SNR teorico
end

figure;
hold on;
plot(bit_values, SNR_empirico, 'o-', 'DisplayName', 'SNR Experimental');
plot(bit_values, SNR_teorico, 'x--', 'DisplayName', 'SNR Teórico');
xlabel('Número de Bits');
ylabel('SNR (dB)');
title('Comparación de SNR Experimental y Teórico en DPCM');
legend;
grid on;
hold off;

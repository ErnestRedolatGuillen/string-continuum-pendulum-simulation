%% Caso continuo (ángulos pequeños)
    % =====================================================================
    % PARAMETROS FÍSICOS
    % =====================================================================
    clear vars;
    clc;
    close all;
    g_grav = 9.81;  % Aceleración de la gravedad (m/s^2)
    L = 1.0;        % Longitud de la cuerda (m)
    
    % =====================================================================
    % CONFIGURACIÓN DE LAS CONDICIONES INICIALES
    % =====================================================================
    % f(y): Posición inicial. Ejemplo: Cuerda desplazada en el centro (tipo pulso)
    f_inicial = @(y) 0.01 * exp(-15 * (y - L/2).^2); 
    
    % g(y): Velocidad inicial. Ejemplo: Perfil de velocidad tras un impacto: -0.2 * (y / L);
    v_inicial = @(y) zeros(size(y));  
    
    % =====================================================================
    % PARÁMETROS NUMÉRICOS
    % =====================================================================
    N_modos = 25;       % Mayor número de modos para capturar transitorios complejos
    N_cc = 64;          % Número de intervalos de Chebyshev (debe ser par)
    N_visual = 200;     % Puntos para la renderización de la cuerda
    t_max = 20;          % Tiempo total de simulación (s)
    dt = 0.05;         % Paso de tiempo de la animación (s)
    
    y_vis = linspace(0, L, N_visual);
    t = 0:dt:t_max;
    
    % =====================================================================
    % 1. CÁLCULO DE RAÍCES DE BESSEL (alpha_n) Y FRECUENCIAS (omega_n)
    % =====================================================================
    alpha = zeros(N_modos, 1);
    omega = zeros(N_modos, 1);
    for n = 1:N_modos
        estimacion = pi * (n - 0.25);
        alpha(n) = fzero(@(z) besselj(0, z), estimacion);
        omega(n) = (alpha(n) / 2) * sqrt(g_grav / L);
    end
    
    % =====================================================================
    % 2. OBTENCIÓN DE NODOS Y PESOS DE CLENSHAW-CURTIS
    % =====================================================================
    [y_cc, w_cc] = pesos_clenshaw_curtis(N_cc, L);
    
    % =====================================================================
    % 3. CÁLCULO DE COEFICIENTES An Y Bn USANDO CUADRATURA
    % =====================================================================
    A = zeros(N_modos, 1);
    B = zeros(N_modos, 1);
    
    % Evaluamos las funciones iniciales en los nodos de Chebyshev
    f_nodos = f_inicial(y_cc);
    v_nodos = v_inicial(y_cc);
    
    % Asegurar que sean vectores columna para evitar conflictos de dimensiones
    f_nodos = f_nodos(:);
    v_nodos = v_nodos(:);
    w_cc = w_cc(:);
    y_cc = y_cc(:);
    
    for n = 1:N_modos
        % Kernel espacial de Bessel para el modo n
        J0_kernel = besselj(0, alpha(n) * sqrt((L - y_cc) / L));
        
        % Denominador común por ortogonalidad
        denominador = L * (besselj(1, alpha(n)))^2;
        
        % Integración por Clenshaw-Curtis (Suma pesada / producto escalar)
        integral_A = sum(w_cc .* f_nodos .* J0_kernel);
        integral_B = sum(w_cc .* v_nodos .* J0_kernel);
        
        % Asignación de coeficientes espectrales
        A(n) = integral_A / denominador;
        B(n) = integral_B / (omega(n) * denominador);
    end
    
    % =====================================================================
    % VÍDEO
    % =====================================================================
    opcion_video = questdlg('¿Deseas guardar la animación en un archivo de video MP4?', ...
	'Guardar Video', ...
	'Sí','No','No'); % 'No' es la opción por defecto

    guardar_video = strcmp(opcion_video, 'Sí');
    
    if guardar_video
    nombre_archivo = sprintf('Sistema_Continuo.mp4');
    v = VideoWriter(nombre_archivo, 'MPEG-4');
    v.FrameRate = 1/dt; % fotogramas por segundo
    v.Quality = 95;   % Alta calidad de compresión
    open(v);
    disp(['Grabando video... El archivo se guardará como: ', nombre_archivo]);
    end
    
    % =====================================================================
    % 4. EVALUACIÓN DE LA SOLUCIÓN Y BUCLE DE ANIMACIÓN
    % =====================================================================
    figure('Color', [1 1 1], 'Position', [200, 100, 550, 650]);
    
    % Precalculamos la matriz espacial de modos para acelerar el bucle
    Modos_Matriz = zeros(N_visual, N_modos);
    for n = 1:N_modos
        Modos_Matriz(:, n) = besselj(0, alpha(n) * sqrt((L - y_vis) / L));
    end
    
    for k = 1:length(t)
        tk = t(k);
        tic;
        % Combinación lineal de las soluciones temporales
        temp_A = A .* cos(omega * tk);
        temp_B = B .* sin(omega * tk);
        componente_temporal = temp_A + temp_B;
        
        % Multiplicación matricial para obtener x(y) en el instante t_k
        x_yt = Modos_Matriz * componente_temporal;
        
        % Gráfica
        plot(x_yt, -y_vis, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2.5);
        hold on;
        plot(0, 0, 'ks', 'MarkerFaceColor', [0.6350, 0.0780, 0.1840], 'MarkerSize', 8); % Soporte
        hold off;
        
        % Estética del gráfico
        grid on;
        xlim([-0.25, 0.25]);
        ylim([-L - 0.1, 0.1]);
        xlabel('Desplazamiento Horizontal x [m]');
        ylabel('Longitud de la Cuerda y [m] (0 = Techo)');
        title(sprintf('Cuerda Continua\nt = %.2f s', tk));
        
        drawnow limitrate;

        if guardar_video
        frame = getframe(gcf);
        writeVideo(v, frame);
        end

        t_p=toc;
        pause(dt-t_p);
    end

if guardar_video
    close(v);
    disp('¡Video guardado y exportado con éxito!');
end

% =========================================================================
% FUNCIÓN AUXILIAR: Cuadratura de Clenshaw-Curtis
% =========================================================================
function [y_cc, w_cc] = pesos_clenshaw_curtis(N, L)
    % Forzar a que N sea par
    if mod(N, 2) ~= 0, N = N + 1; end
    
    % 1. Nodos de Chebyshev en [-1, 1]
    theta = pi * (0:N)' / N;
    x_cheb = cos(theta);
    
    % Mapeo de nodos al intervalo físico de la cuerda [0, L]
    y_cc = (L / 2) * (x_cheb + 1);
    
    % 2. Cálculo directo de los pesos estándar en [-1, 1]
    w_std = zeros(N + 1, 1);
    for i = 0:N
        sum_terms = 0;
        for j = 1:(N/2 - 1)
            sum_terms = sum_terms + (2 / (1 - 4*j^2)) * cos(2*j*i*pi / N);
        end
        % Añadir los términos extremos de la sumatoria analítica
        val_izq = 1;
        val_der = cos(i*pi) / (1 - N^2);
        
        w_std(i+1) = (2 / N) * (0.5 * val_izq + sum_terms + 0.5 * val_der);
    end
    
    % Ajustar extremos del intervalo
    w_std(1) = w_std(1) / 2;
    w_std(end) = w_std(end) / 2;
    
    % Escalar pesos al tamaño del dominio [0, L]
    w_cc = w_std * (L / 2);
    
    % Invertir vectores para que comiencen físicamente en y = 0 (anclaje fijo)
    y_cc = flipud(y_cc);
    w_cc = flipud(w_cc);
end
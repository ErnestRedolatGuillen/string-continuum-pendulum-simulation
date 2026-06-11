%% Animación Simultánea de Modos Normales (del sistema lineal)
clear; clc; close all;

% PARÁMETROS FÍSICOS
n = 6;              % # péndulos
L = 1;              % Longitud total (m)
M_total = 20;       % Masa total (kg)
g = 9.81;           % Gravedad (m/s^2)
l = L/n;            % Longitud de cada segmento
m = M_total/n;      % Masa de cada partícula

% 1. CONSTRUCCIÓN DE MATRICES (M y K)
K = zeros(n,n);
M = zeros(n,n);
for i = 1:n
    K(i,i) = m * g * l * (n - i + 1); 
    for j = 1:n
        M(i,j) = m * l^2 * (n - max(i,j) + 1);
    end
end

% 2. CÁLCULO DE MODOS PROPIOS (VAPs y VEPs)
[V, D] = eig(K, M);
frecuencias_angulares = sqrt(diag(D));

% Ordenar modos de menor a mayor frecuencia
[frecuencias_angulares, idx] = sort(frecuencias_angulares);
V = V(:, idx);

% =========================================================================
% PREGUNTA INTERACTIVA PARA GUARDAR VIDEO
% =========================================================================
opcion_video = questdlg('¿Deseas guardar la animación en un archivo de video MP4?', ...
	'Guardar Video', ...
	'Sí','No','No'); % 'No' es la opción por defecto

guardar_video = strcmp(opcion_video, 'Sí');

if guardar_video
    nombre_archivo = sprintf('Modos_Normales_n%d.mp4', n);
    v = VideoWriter(nombre_archivo, 'MPEG-4');
    v.FrameRate = 50; % 50 fotogramas por segundo
    v.Quality = 95;   % Alta calidad de compresión
    open(v);
    disp(['Grabando video... El archivo se guardará como: ', nombre_archivo]);
end
% =========================================================================

% 3. CONFIGURACIÓN DE LA FIGURA Y SUBPLOTS
fig = figure('Color', 'w', 'Position', [100, 100, 900, 600]);

modos_a_graficar = min(n, 12); 
columnas = ceil(sqrt(modos_a_graficar));
filas = ceil(modos_a_graficar / columnas);

amplitud_visual = 0.25; 
hCuerdas = cell(modos_a_graficar, 1); 

% Dibujado inicial de los subplots
for modo = 1:modos_a_graficar
    subplot(filas, columnas, modo);
    grid on; hold on; axis equal;
    
    % Línea de equilibrio de fondo
    plot([0, 0], [0, -L], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    
    % Crear el objeto gráfico animable
    hCuerdas{modo} = plot(zeros(n+1,1), -(0:l:L), '-or', ...
        'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    
    xlim([-L*0.6, L*0.6]);
    ylim([-L*1.1, L*0.1]);
    
    f_hz = frecuencias_angulares(modo) / (2*pi);
    title(sprintf('Modo %d (%.2f s-1)', modo, f_hz), 'FontSize', 10);
    
    if modo > (filas-1)*columnas, xlabel('X (m)'); end
    if mod(modo, columnas) == 1, ylabel('Y (m)'); end
end

sgtitle(sprintf('Movimiento de los Modos Normales (n = %d)', n), ...
    'FontSize', 14, 'FontWeight', 'bold');

% 4. BUCLE DE ANIMACIÓN
t_anim = 0:0.02:10;
fps_control = tic;

for k = 1:length(t_anim)
    t = t_anim(k);
    
    % Actualizar cada modo individualmente
    for modo = 1:modos_a_graficar
        theta_actual = amplitud_visual * (V(:, modo) / max(abs(V(:, modo)))) * cos(frecuencias_angulares(modo) * t);
        
        posX = [0; cumsum(l * sin(theta_actual))];
        posY = [0; cumsum(-l * cos(theta_actual))];
        
        set(hCuerdas{modo}, 'XData', posX, 'YData', posY);
    end
    
    % Forzar el renderizado en pantalla
    drawnow;
    
    % Captura el frame actual de la figura
    if guardar_video
        frame = getframe(fig);
        writeVideo(v, frame);
    end
    
    % Sincronización para la reproducción en tiempo real
    t_procesado = toc(fps_control);
    pause(max(0, 0.02 - t_procesado));
    fps_control = tic; 
end

% Cierre del archivo de video
if guardar_video
    close(v);
    disp('¡Video guardado y exportado con éxito!');
end
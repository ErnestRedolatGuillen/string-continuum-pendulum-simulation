%% Caso discreto para ángulos pequeños
close all;
% PARÁMETROS FÍSICOS
n = 20;             % Número de péndulos
L = 1;              % Longitud total (m)
M_total = 20;      % Masa total (kg)
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

% 2. CONFIGURACIÓN DE LA SIMULACIÓN
tspan = 0:0.05:20;  % Tiempo de 0 a 30s con pasos de 0.05s
theta0 = [0.17*ones(n,1)]; % Posiciones iniciales
omega0 = [zeros(n,1)]; %Velocidades iniciales

y0 = [theta0; omega0];    % Estado inicial [ángulos; velocidades]

% Definición del sistema de EDOs: y' = [velocidades; aceleraciones]
sistema = @(t, y) [y(n+1:end); -M \ (K * y(1:n))];

% Resolver las ecuaciones
[t, y_out] = ode45(sistema, tspan, y0);
thetas = y_out(:, 1:n);


%% VIDEO
opcion_video = questdlg('¿Deseas guardar la animación en un archivo de video MP4?', ...
	'Guardar Video', ...
	'Sí','No','No'); % 'No' es la opción por defecto

guardar_video = strcmp(opcion_video, 'Sí');

if guardar_video
    nombre_archivo = sprintf('Sistema_Lineal_n%d.mp4', n);
    v = VideoWriter(nombre_archivo, 'MPEG-4');
    v.FrameRate = (size(tspan,2)-1)/tspan(end); % fotogramas por segundo
    v.Quality = 95;   % Alta calidad de compresión
    open(v);
    disp(['Grabando video... El archivo se guardará como: ', nombre_archivo]);
end
%%

% 3. ANIMACIÓN
figure('Color', 'w');
grid on;
hold on;
axis equal;
% Ajustar los límites del gráfico
xlim([-L*0.6, L*0.6]); 
ylim([-L*1.1, L*0.1]);
xlabel('X (m)'); ylabel('Y (m)');


% Crear el objeto gráfico de la cuerda (línea con círculos)
hCuerda = plot(0, 0, '-or', 'LineWidth', 2, 'MarkerFaceColor', 'r');

% Bucle de frames para el movimiento
for k = 1:length(t)
    tk=tspan(k);
    title(sprintf('Simulación del sistema lineal (n =%5.0f)\nt = %5.2f s',n,tk));
    tic;
    % Calcular posiciones cartesianas (x, y) a partir de los ángulos
    % x = l * sin(theta1) + l * sin(theta2) ...
    current_thetas = thetas(k, :);
    posX = [0, cumsum(l * sin(current_thetas))];
    posY = [0, cumsum(-l * cos(current_thetas))];
    
    % Actualizar datos del gráfico
    set(hCuerda, 'XData', posX, 'YData', posY);
    
    % Forzar dibujado y pequeña pausa (que depende del tiempo de procesado) para velocidad real
    drawnow;

    % Captura el frame actual de la figura
    if guardar_video
        frame = getframe(gcf);
        writeVideo(v, frame);
    end
    tproces=toc;
    pause(max(0,0.05-tproces)); 
end

if guardar_video
    close(v);
    disp('¡Video guardado y exportado con éxito!');
end
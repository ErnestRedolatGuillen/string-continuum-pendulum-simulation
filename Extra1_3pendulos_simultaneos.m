%% Caso_discreto_triple
% Simulación simultánea (en la misma gráfica) de tres sistemas de péndulos
% múltiples con condiciones iniciales personalizables para cada sistema.
% Útil para mostrar la sensibilidad a las CI
close all;
% --- Configuración del Sistema ---
n = 20;             % Número de péndulos por cadena
L_total = 1.0;      % Longitud total (metros)
M_total = 20.0;     % Masa total (kg)
g = 9.81;           

l = L_total / n;
m = M_total / n;

% --- Condiciones Iniciales---
dx=1e-8;

% Péndulo 1:
theta0_1 = deg2rad(90+dx) * ones(n, 1); 

% Péndulo 2: 
theta0_2 = deg2rad(90) * ones(n, 1);%

% Péndulo 3:
theta0_3 = deg2rad(90-dx) * ones(n, 1);

% Velocidades iniciales (todas a cero)
omega0 = zeros(n, 1); 

% Vector de estado global: [angulos1; angulos2; angulos3; veloc1; veloc2; veloc3]
z0 = [theta0_1; theta0_2; theta0_3; omega0; omega0; omega0];

% --- Simulación ---
tspan = 0:0.05:20; % Pasos de tiempo fijos para la animación
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Llamada al solver con el sistema adaptado para 3 péndulos
[t, z] = ode45(@(t, z) sistema_derivadas_triple(t, z, n, l, m, g), tspan, z0, opts);

% VÍDEO
% 1. Preparar el archivo de video
opcion_video = questdlg('¿Deseas guardar la animación en un archivo de video MP4?', ...
	'Guardar Video', ...
	'Sí','No','No'); % 'No' es la opción por defecto

guardar_video = strcmp(opcion_video, 'Sí');

if guardar_video
v = VideoWriter('discreto_triple.mp4', 'MPEG-4');
v.FrameRate = (size(tspan,2)-1)/tspan(end); 
open(v);
disp(['Grabando video... El archivo se guardará como: ', nombre_archivo]);
end

% --- Animación ---
figure('Color', 'w', 'Name', sprintf('Simulación Simultánea de 3 Péndulos (n=%d)', n));
hold on;
% Creamos tres plots con colores diferentes
h_plot1 = plot(0, 0, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'Color', 'r', 'DisplayName', 'Péndulo 1');
h_plot2 = plot(0, 0, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'Color', 'b', 'DisplayName', 'Péndulo 2');
h_plot3 = plot(0, 0, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'g', 'Color', 'g', 'DisplayName', 'Péndulo 3');
grid on;
axis equal;
legend('Location', 'northeast');

% Ajustar ejes según longitud total
limit = L_total * 1.1;
axis([-limit/2, limit/2, -limit, 0.2]);
xlabel('X (m)'); ylabel('Y (m)');


for k = 1:length(t)
    tic;
    
    % Extraer ángulos de cada péndulo en el instante k
    thetas_1 = z(k, 1:n)';
    thetas_2 = z(k, n+1:2*n)';
    thetas_3 = z(k, 2*n+1:3*n)';
    
    % Calcular posiciones cartesianas para el Péndulo 1
    x1 = [0; cumsum(l * sin(thetas_1))];
    y1 = [0; cumsum(-l * cos(thetas_1))];
    
    % Calcular posiciones cartesianas para el Péndulo 2
    x2 = [0; cumsum(l * sin(thetas_2))];
    y2 = [0; cumsum(-l * cos(thetas_2))];
    
    % Calcular posiciones cartesianas para el Péndulo 3
    x3 = [0; cumsum(l * sin(thetas_3))];
    y3 = [0; cumsum(-l * cos(thetas_3))];
    
    % Actualizar los tres gráficos simultáneamente
    set(h_plot1, 'XData', x1, 'YData', y1);
    set(h_plot2, 'XData', x2, 'YData', y2);
    set(h_plot3, 'XData', x3, 'YData', y3);
    
    title(sprintf('Péndulos múltiples simultáneos\nTiempo: %.2f s', t(k)));
    drawnow; 
    
    tproces = toc;
    pause(max(0, 0.05 - tproces));
    
    % VÍDEO
    if guardar_video
    % 2. Capturar el cuadro actual de la figura
    frame = getframe(gcf); 
    % 3. Escribir el cuadro en el video
    writeVideo(v, frame);
    end
end

% 4. Cerrar el archivo de vídeo
if guardar_video
close(v);
disp('Video guardado con éxito.');
end

%% Definimos el sistema para 3 péndulos independientes
function dzdt = sistema_derivadas_triple(~, z, n, l, m, g)
    % Desempaquetar ángulos y velocidades para los 3 sistemas
    theta1 = z(1:n);
    theta2 = z(n+1:2*n);
    theta3 = z(2*n+1:3*n);
    
    omega1 = z(3*n+1:4*n);
    omega2 = z(4*n+1:5*n);
    omega3 = z(5*n+1:end);
    
    % Como los pèndulos no interactúan entre sí físicamente (solo se simulan a la vez),
    % calculamos las aceleraciones (alpha) de cada uno por separado usando la misma lógica.
    alpha1 = resolver_pendulo(theta1, omega1, n, l, m, g);
    alpha2 = resolver_pendulo(theta2, omega2, n, l, m, g);
    alpha3 = resolver_pendulo(theta3, omega3, n, l, m, g);
    
    % Empaquetar de nuevo el vector de derivadas: [velocidades; aceleraciones]
    dzdt = [omega1; omega2; omega3; alpha1; alpha2; alpha3];
end

% Función auxiliar para evitar duplicar el código de Euler-Lagrange
function alpha = resolver_pendulo(theta, omega, n, l, m, g)
    M_mat = zeros(n, n);
    F = zeros(n, 1);
    
    for m_idx = 1:n
        F(m_idx) = -m * g * l * (n - m_idx + 1) * sin(theta(m_idx));
        
        for j = 1:n
            A_mj = (n - max(m_idx, j) + 1) * m * l^2;
            M_mat(m_idx, j) = A_mj * cos(theta(m_idx) - theta(j));
            
            if m_idx ~= j
                F(m_idx) = F(m_idx) - A_mj * omega(j)^2 * sin(theta(m_idx) - theta(j));
            end
        end
    end
    alpha = M_mat \ F;
end
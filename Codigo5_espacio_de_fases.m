function Codigo5_espacio_de_fases()
% Representación del espacio de fases 3D del doble péndulo considerando la 
% energía constante
close all;
    % --- PARÁMETROS FÍSICOS ---
    g = 9.81;
    l = 1.0;
    m = 20.0;
    
    % --- CONDICIONES INICIALES ---
    theta1_0 = deg2rad(80);
    theta2_0 = deg2rad(0);
    omega1_0 = 0;
    omega2_0 = 0;
    z = [theta1_0; theta2_0; omega1_0; omega2_0]; % Vector de estado inicial
    
    % --- CONFIGURACIÓN DE LA ANIMACIÓN ---
    t_total = 150;        % Tiempo total de simulación
    dt = 0.05;           % Intervalo de tiempo para cada paso de animación

    % ---Vídeo
    opcion_video = questdlg('¿Deseas guardar la animación en un archivo de video MP4?', ...
	'Guardar Video', ...
	'Sí','No','No'); % 'No' es la opción por defecto

    guardar_video = strcmp(opcion_video, 'Sí');

    if guardar_video
    nombre_archivo = sprintf('Espacio_Fases.mp4');
    v = VideoWriter(nombre_archivo, 'MPEG-4');
    v.FrameRate = 1/dt; % fotogramas por segundo
    v.Quality = 95;   % Alta calidad de compresión
    open(v);
    disp(['Grabando video... El archivo se guardará como: ', nombre_archivo]);
    end
    % ---
    
    % Configuración de la ventana gráfica
    figure('Position', [200, 100, 800, 600], 'Color', 'w');
    
    % Creamos una línea animada en 3D (el camino del espacio de fases)
    h_trayectoria = animatedline('Color', [0 0.5 0.8], 'LineWidth', 1.5);
    
    % Creamos un punto que representa el "estado actual" del péndulo
    hold on;
    h_estado_actual = plot3(z(1), z(2), z(4), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    hold off;
    
    % Configuración de los ejes y etiquetas
    xlabel('\theta_1 (Posición Péndulo 1) [rad]');
    ylabel('\theta_2 (Posición Péndulo 2) [rad]');
    zlabel('\omega_2 (Velocidad Péndulo 2) [rad/s]');
    title('Espacio de Fases 3D en Tiempo Real');
    
    % Límites fijos para evitar que los ejes salten durante la animación
    xlim([-pi, pi]);
    ylim([-pi, pi]);
    zlim([-10, 10]); 
    
    grid on;
    view(45, 30); % Cambiamos la perspectiva de la cámara 3D
    rotate3d on;  % Permite hacer clic y rotar el gráfico mientras corre
    
    % Opciones de precisión para el integrador
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    
    % --- BUCLE DE SIMULACIÓN EN DIRECTO ---
    fprintf('Iniciando animación... Puedes rotar el gráfico con el ratón.\n');
    
    for t_actual = 0:dt:t_total
        % Si la ventana gráfica se cierra, detenemos el programa
        if ~ishandle(h_trayectoria)
            break;
        end
        
        % Resolver el sistema solo para el pequeño intervalo [t_actual, t_actual + dt]
        [~, z_sol] = ode45(@(t, z) derivadas(t, z, g, l, m), [t_actual, t_actual + dt], z, options);
        
        % El último punto de esta solución será el estado actual para la siguiente iteración
        z = z_sol(end, :)';
        
        % Normalizar los ángulos al rango [-pi, pi] para que no se escapen del gráfico
        theta1_norm = mod(z(1) + pi, 2*pi) - pi;
        theta2_norm = mod(z(2) + pi, 2*pi) - pi;
        omega2_actual = z(4);
        
        % Añadir el nuevo punto a la línea animada
        addpoints(h_trayectoria, theta1_norm, theta2_norm, omega2_actual);
        
        % Actualizar la posición del punto rojo indicador
        set(h_estado_actual, 'XData', theta1_norm, 'YData', theta2_norm, 'ZData', omega2_actual);
        
        % Forzar a Matlab a dibujar los cambios en la pantalla
        drawnow;

        if guardar_video
        frame = getframe(gcf);
        writeVideo(v, frame);
        end
    end
    fprintf('Animación finalizada.\n');

if guardar_video
    close(v);
    disp('¡Video guardado y exportado con éxito!');
end
end

% --- MOVIMIENTO DEL PÉNDULO DOBLE ---
function dzdt = derivadas(~, z, g, l, m)
    t1 = z(1); t2 = z(2);
    w1 = z(3); w2 = z(4);
    
    % Matriz de inercia
    M = [2*m*l^2,          m*l^2*cos(t1 - t2);
         m*l^2*cos(t1 - t2), m*l^2];
         
    % Fuerzas centrífugas y gravedad
    F = [-2*m*g*l*sin(t1) - m*l^2*w2^2*sin(t1 - t2);
         -m*g*l*sin(t2)   + m*l^2*w1^2*sin(t1 - t2)];
         
    % Aceleraciones
    alpha = M \ F;
    
    dzdt = [w1; w2; alpha(1); alpha(2)];
end
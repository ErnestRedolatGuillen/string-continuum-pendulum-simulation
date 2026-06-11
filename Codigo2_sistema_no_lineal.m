%% Caso_discreto_no_lineal
close all;
    % --- Configuración del Sistema ---
    n = 20;             % Número de péndulos
    L_total = 1.0;      % Longitud total (metros)
    M_total = 20.0;      % Masa total (kg)
    g = 9.81;           
    
    l = L_total / n;
    m = M_total / n;
    
    % --- Condiciones Iniciales ---

    % Ángulos iniciales: theta_i_0 en radianes
    % Ejemplos: una curva suave inicial (senoidal) 0.5 * sin(linspace(0, pi/2, n))'; 
    %theta0=zeros(n,1); theta0(1) = deg2rad(45); % Solo el primero a 45 grados

    theta0 = deg2rad(90)*ones(n, 1);% Todos en 90º
    omega0 = zeros(n, 1); % Todas en 0 s-1
    z0 = [theta0; omega0];
    
    % --- Simulación ---
    tspan = 0:0.05:20; % Pasos de tiempo fijos para la animación
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    % El sistema está programado en una función aparte.
    [t, z] = ode45(@(t, z) sistema_derivadas(t, z, n, l, m, g), tspan, z0, opts);
    
    %--- Vídeo ---
    opcion_video = questdlg('¿Deseas guardar la animación en un archivo de video MP4?', ...
	'Guardar Video', ...
	'Sí','No','No'); % 'No' es la opción por defecto

    guardar_video = strcmp(opcion_video, 'Sí');
    
    if guardar_video
    nombre_archivo = sprintf('Sistema__No_Lineal_n%d.mp4', n);
    v = VideoWriter(nombre_archivo, 'MPEG-4');
    v.FrameRate = (size(tspan,2)-1)/tspan(end); % fotogramas por segundo
    v.Quality = 95;   % Alta calidad de compresión
    open(v);
    disp(['Grabando video... El archivo se guardará como: ', nombre_archivo]);
    end

    % --- Animación ---
    figure('Color', 'w', 'Name', sprintf('Simulación n=%5.0f',n));
    h_plot = plot(0, 0, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    grid on;
    axis equal;
    % Ajustar ejes según longitud total
    limit = L_total * 1.1;
    axis([-limit/2, limit/2, -limit, 0.2]);
    xlabel('X (m)'); ylabel('Y (m)');
    
    for k = 1:length(t)
        tic;
        % Extraer ángulos en el instante k
        thetas = z(k, 1:n);
        
        % Calcular posiciones cartesianas de cada masa (sumatoria acumulada)
        x = [0; cumsum(l * sin(thetas'))];
        y = [0; cumsum(-l * cos(thetas'))];
        
        % Actualizar gráfico
        set(h_plot, 'XData', x, 'YData', y);
        title(sprintf('Simulación del sistema no lineal\nTiempo: %.2f s', t(k)));
        drawnow;

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

%% Definimos el sistema en una función aparte
function dzdt = sistema_derivadas(~, z, n, l, m, g)
    theta = z(1:n);
    omega = z(n+1:end);
    
    M_mat = zeros(n, n);
    F = zeros(n, 1);
    
    % Llenado de matrices basado en Euler-Lagrange
    for m_idx = 1:n
        % Fuerza restauradora de la gravedad
        F(m_idx) = -m * g * l * (n - m_idx + 1) * sin(theta(m_idx));
        
        for j = 1:n
            % Coeficiente de inercia acoplada
            % Representa la masa total que cuelga por debajo de los nodos m y j
            A_mj = (n - max(m_idx, j) + 1) * m * l^2;
            
            % Matriz de Masa: M(theta) * ddot{theta}
            M_mat(m_idx, j) = A_mj * cos(theta(m_idx) - theta(j));
            
            % Términos centrífugos: A_mj * dot{theta}^2 * sin(theta_m - theta_j)
            if m_idx ~= j
                F(m_idx) = F(m_idx) - A_mj * omega(j)^2 * sin(theta(m_idx) - theta(j));
            end
        end
    end
    
    % Resolver sistema lineal para aceleraciones
    alpha = M_mat \ F;
    dzdt = [omega; alpha];
end
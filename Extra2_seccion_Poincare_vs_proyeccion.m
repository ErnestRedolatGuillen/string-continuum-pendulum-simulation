function Extra2_seccion_Poincare_vs_proyeccion()
% Comparación de la visualización de la proyección del espacio de fases 3D
% de un péndulo doble en el plano frente a su sección de Poincaré, con CI
% configurable.

% Dependiendo del ángulo (energía), la sección de Poincaré muestra patrones
% ordenados o dispersos.

    % --- PARÁMETROS FÍSICOS --- 
    g = 9.81;
    l = 1;
    m = 20.0;
    
    % --- CONDICIONES INICIALES --- (doble péndulo)
    theta1_0 = deg2rad(50); % Primer péndulo a 90 grados
    theta2_0 = deg2rad(0);  % Segundo péndulo vertical
    omega1_0 = 0;
    omega2_0 = 0;
    z0 = [theta1_0; theta2_0; omega1_0; omega2_0];
    
    % --- TIEMPO DE SIMULACIÓN (Largo para acumular puntos en el fractal) ---
    tspan = 0:0.05:2000;
    
    % Opciones de precisión alta para evitar que el error numérico destruya el fractal
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    
    fprintf('Simulando el sistema... Por favor espera.\n');
    [t, z] = ode45(@(t, z) derivadas(t, z, g, l, m), tspan, z0, options);
    fprintf('Simulación completada. Procesando la Sección de Poincaré...\n');
    
    % --- EXTRACCIÓN DE LA SECCIÓN DE POINCARÉ ---
    % Cortamos el espacio de fases cuando el péndulo 1 pasa por la vertical (theta1 = 0)
    % con velocidad positiva (omega1 > 0) para evitar duplicidad de espejo.
    
    theta1 = z(:, 1);
    omega1 = z(:, 3);
    theta2 = z(:, 2);
    omega2 = z(:, 4);
    
    % Encontrar los índices donde theta1 cruza por cero
    cruces_indices = [];
    for i = 1:(length(theta1) - 1)
        if theta1(i) * theta1(i+1) < 0 && omega1(i) > 0
            % Interpolación lineal simple para mayor precisión en el punto de cruce
            frac = -theta1(i) / (theta1(i+1) - theta1(i));
            t_cruce = t(i) + frac * (t(i+1) - t(i));
            theta2_cruce = theta2(i) + frac * (theta2(i+1) - theta2(i));
            omega2_cruce = omega2(i) + frac * (omega2(i+1) - omega2(i));
            
            % Normalizar theta2 al rango [-pi, pi] debido a las rotaciones de 360°
            theta2_cruce = mod(theta2_cruce + pi, 2*pi) - pi;
            
            cruces_indices = [cruces_indices; theta2_cruce, omega2_cruce];
        end
    end
    
    % --- GRAFICAR LOS RESULTADOS ---
    figure('Position', [100, 100, 1100, 500]);
    
    % Gráfico 1: Espacio de Fases Continuo proyectado sobre el plano theta1=0
    subplot(1, 2, 1);
    plot(mod(theta2 + pi, 2*pi) - pi, omega2, 'Color', [0 0.4470 0.7410 0.1], 'LineWidth', 0.5);
    xlabel('\theta_2 (rad)');
    ylabel('\omega_2 (rad/s)');
    title('Espacio de Fases Continuo proyectado en el plano theta1=0 (Péndulo 2)');
    grid on;
    
    % Gráfico 2: Sección de Poincaré (intersección de las órbitas con el
    % plano theta1=0 y cuando tienen omega1>0
    subplot(1, 2, 2);
    if ~isempty(cruces_indices)
        plot(cruces_indices(:, 1), cruces_indices(:, 2), '.k', 'MarkerSize', 3);
    else
        text(0.5, 0.5, 'No se detectaron cruces', 'HorizontalAlignment', 'center');
    end
    xlabel('\theta_2 (rad)');
    ylabel('\omega_2 (rad/s)');
    title('Sección de Poincaré (\theta_1 = 0, \omega_1 > 0)');
    grid on;
    
end

% --- FUNCIÓN DE DERIVADAS NO LINEALES (Formulada mediante Lagrange) ---
function dzdt = derivadas(~, z, g, l, m)
    % Desempaquetar variables de estado
    t1 = z(1);
    t2 = z(2);
    w1 = z(3);
    w2 = z(4);
    
    % Matriz de Masa M(theta) original no lineal para n=2
    M = [2*m*l^2,          m*l^2*cos(t1 - t2);
         m*l^2*cos(t1 - t2), m*l^2];
         
    % Vector de Fuerzas F(theta, omega) que incluye gravedad y fuerzas centrífugas
    F = [-2*m*g*l*sin(t1) - m*l^2*w2^2*sin(t1 - t2);
         -m*g*l*sin(t2)   + m*l^2*w1^2*sin(t1 - t2)];
         
    % Resolver para las aceleraciones angulares (alpha = M \ F)
    alpha = M \ F;
    
    % Vector de derivadas de primer orden: [velocidades; aceleraciones]
    dzdt = [w1; w2; alpha(1); alpha(2)];
end
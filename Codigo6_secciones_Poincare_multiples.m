function Codigo6_secciones_Poincare_multiples()
% Representación de la secciones de Poincaré múltiples para el péndulo
% doble. En cada subplot se muestran varias secciones de Poincaré
% superpuestas de distintas trayectorias (cada una de un mismo color) que
% tienen la misma energía.
close all;
    % --- Configuración de parámetros físicos (n=2) ---
    g = 9.81; 
    l1 = 0.5; l2 = 0.5;
    m1 = 10.0; m2 = 10.0;
    
    % --- REFERENCIA DE ENERGÍA ---
    % E = 0 J  en reposo.
    E_deseada = [2 15 40 70 100 150];
    
    num_trayectorias = 50; 
    tspan = [0, 150]; 
    
    figure('Color', 'w', 'Position', [100, 100, 800, 600]);
    hold on; grid on; box on;
    
    fprintf('Calculando Sección de Poincaré para Energía = %.2f J', E_deseada);
    
    % Rango de muestreo para theta1 (ajustado para energías bajas)
    theta1_vals = linspace(-0.4, 0.4, num_trayectorias); 

    sk=0;
    sgtitle('Secciones de Poincaré','FontSize',20);

    for E=E_deseada
    sk=sk+1;
    subplot(2,3,sk);
    hold on;
    title(sprintf('Energía = %.0f J', E), 'FontSize', 14);
    xlabel('\theta_1 (rad)', 'FontSize', 12);
    ylabel('\omega_1 (rad/s)', 'FontSize', 12);


    for i = 1:num_trayectorias
        th1 = theta1_vals(i);
        th2 = 0;
        om2 = 0;
        
        % Energía Potencial MODIFICADA: Se le suma el desfase para que en (0,0) sea V = 0
        V_original = -m1*g*l1*cos(th1) - m2*g*(l1*cos(th1) + l2*cos(th2));
        V_offset = (m1 + m2)*g*l1 + m2*g*l2; % Energía potencial en reposo
        V = V_original + V_offset; 
        
        % Energía Cinética requerida
        T_req = E - V;
        
        if T_req < 0
            continue; % Energía inalcanzable para este ángulo inicial
        end
        
        % Velocidad angular inicial om1
        om1 = sqrt(2 * T_req / ((m1 + m2) * l1^2));
        z0 = [th1; th2; om1; om2];
        
        % Configuración del detector de eventos (Sección de Poincaré)
        opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9, 'Events', @event_poincare);
        
        % Simulación
        [~, ~, te, ye, ie] = ode45(@(t, z) ec_movimiento(z, m1, m2, l1, l2, g), tspan, z0, opts);
        
        % Filtrar y graficar
        if ~isempty(ye)
            validos = ye(:,4) > 0; % Condición omega2 > 0
            th1_puntos = ye(validos, 1);
            om1_puntos = ye(validos, 3);
            
            % Normalizar ángulo
            th1_puntos = mod(th1_puntos + pi, 2*pi) - pi;
            plot(th1_puntos, om1_puntos, '.', 'MarkerSize', 4);
        end
    end

    end
    
end

%% --- FUNCIÓN DE EVENTOS ---
function [value, isterminal, direction] = event_poincare(~, z)
    value = z(2); % Detectar theta2 = 0
    isterminal = 0; 
    direction = 0;  
end

%% --- ECUACIONES DE MOVIMIENTO ---
% Nota: Las ecuaciones de movimiento NO cambian, ya que las fuerzas dependen
% de las derivadas del potencial (gradientes), y la constante que sumamos desaparece al derivar.
function dzdt = ec_movimiento(z, m1, m2, l1, l2, g)
    th1 = z(1); th2 = z(2);
    om1 = z(3); om2 = z(4);
    
    M = [ (m1 + m2)*l1^2,          m2*l1*l2*cos(th1 - th2);
          m2*l1*l2*cos(th1 - th2), m2*l2^2 ];
      
    F = [ -m2*l1*l2*om2^2*sin(th1 - th2) - (m1 + m2)*g*l1*sin(th1);
           m2*l1*l2*om1^2*sin(th1 - th2) - m2*g*l2*sin(th2) ];
       
    alpha = M \ F;
    dzdt = [om1; om2; alpha(1); alpha(2)];
end
%% Define All Parameters at the Start
% Ambient Conditions
Ta = 21;           % Ambient temperature (°C)
V = 2;              % Air velocity (m/s)

% Geometric Parameters
L_s = 35e-3;        % Characteristic length for surface (m)
L_c = 1.8;          % Characteristic length for tube (m)
beta_deg = 10;      % Inclination angle (degrees)

% Material Properties
epsilon_s = 0.95;   % Emissivity of the surface
epsilon_c = 0.89;   % Emissivity of the tube

% Physical Constants
g = 9.81;           % Acceleration due to gravity (m/s^2)
sigma = 5.67e-8;    % Stefan-Boltzmann constant (W/m^2·K^4)

% Loop over all combinations of Ts and Tc
for i = 1:length(Ts_range)
    for j = 1:length(Tc_range)
        % Set initial temperature guesses
        Ts = Ts_range(i);  % Initial surface temperature (°C)
        Tc = Tc_range(j);  % Initial tube temperature (°C)
        
        fprintf('\nTesting Ts = %.1f°C, Tc = %.1f°C\n', Ts, Tc);
        
        % Run the existing calculations
        T_out_result = calculate_T_out(T_out, Ts, Tc, Ta, m_dot, Cp, T_in, G, A_s, alpha_s, V, L_c, air_properties, epsilon_c, sigma, error_boundary, L_s, beta_deg, g, epsilon_s);
        Ts_result = inner(Tc, Ts, T_out_result, air_properties, L_s, beta_deg, g, epsilon_s, epsilon_c, sigma, V, L_c, Ta, m_dot, Cp, T_in, h_i, A_i, error_boundary);
        Tc_result = outer(Tc, Ts_result, air_properties, L_s, beta_deg, g, epsilon_s, epsilon_c, sigma, V, L_c, Ta, A_s, error_boundary);
        
        % Store results
        T_out_results(i, j) = T_out_result;
        Ts_results(i, j) = Ts_result;
        Tc_results(i, j) = Tc_result;
        
        % Display final results for this combination
        fprintf('Final Results: T_out = %.6f°C, Ts = %.6f°C, Tc = %.6f°C\n', T_out_result, Ts_result, Tc_result);
    end
end    % Initial tube temperature (°C) 30

% Optionally, display all results in a table
disp('Summary of Results:');
results_table = table(Ts_range', Ts_results, Tc_results, T_out_results, ...
    'VariableNames', {'Initial_Ts', 'Final_Ts', 'Final_Tc', 'Final_T_out'});
disp(results_table);

error_boundary = 0.000000000001;

air_properties = [...
     5,   0.02401, 1.382e-5, 0.7350;
    10,   0.02439, 1.426e-5, 0.7336;
    15,   0.02476, 1.470e-5, 0.7323;
    20,   0.02514, 1.516e-5, 0.7309;
    25,   0.02551, 1.562e-5, 0.7296;
    30,   0.02588, 1.608e-5, 0.7282;
    35,   0.02625, 1.655e-5, 0.7268;
    40,   0.02662, 1.702e-5, 0.7255;
    50,   0.02735, 1.798e-5, 0.7228;
    60,   0.02808, 1.896e-5, 0.7202;
    70,   0.02881, 1.995e-5, 0.7177;
    80,   0.02953, 2.097e-5, 0.7154;
    90,   0.03024, 2.201e-5, 0.7132;
   100,   0.03095, 2.306e-5, 0.7111];

%% Work out T_out from equation 2

function [T_out_new] = calculate_T_out(T_out, Ts, Tc, Ta, m_dot, Cp, T_in, G, A_s, alpha_s, V, L_c, air_properties, epsilon_c, sigma, error_boundary, L_s, beta_deg, g, epsilon_s)
    abs_error = inf; % Initialize error to a large value
    T_out_new = T_out; % Start with initial guess for T_out

    while abs_error > error_boundary
        % Step 1: Compute Convective Heat Transfer Coefficient (Surface to Cover)
        [~, ~, ~, ~, h_s_c] = heatTransferCalc(Ts, Tc, L_s, beta_deg, g, air_properties);

        % Step 2: Compute Radiative Heat Transfer Coefficient (Surface to Cover)
        hr_s_c = radiativeHeatTransfer(Ts, Tc, epsilon_s, epsilon_c, sigma);
    
        % Step 3: Compute Convective Heat Transfer Coefficient (Cover to Ambient)
        [~, ~, h_c_a] = calculateHeatTransfer(V, L_c, air_properties, Tc, Ta);
    
        % Step 4: Compute Radiative Heat Transfer Coefficient (Cover to Ambient)
        hr_c_a = calculateRadiativeHeatTransfer(Ta, Tc, epsilon_c, sigma);
    
        % Step 5: Compute Overall Heat Transfer Coefficient U
        U_new = ((1 / (h_s_c + hr_s_c)) + (1 / (h_c_a + hr_c_a)))^(-1);
    
        % Step 6: Calculate new T_out using Equation 2 (modified)
        tau_c = 0.92;
        T_out_old = T_out_new;
        % Compute the mean temperature for heat loss
        T_mean = (Ts + Tc) / 2;
        % Compute the heat loss term
        heat_loss = U_new * A_s * (T_mean - Ta);
        % Compute the absorbed solar energy
        Q_solar = G * A_s * tau_c * alpha_s;
        % Compute the denominator
        denom = m_dot * Cp;
        % Compute the temperature difference
        T_out_new = T_in + (Q_solar - heat_loss) / denom;
    
        % Step 7: Compute absolute error
        abs_error = abs(T_out_new - T_out_old);
    
        % Step 8: Display progress
        fprintf('Current T_out = %.6f °C, Error = %.6f\n', T_out_new, abs_error);
    end
end

% Initial guess for T_out
T_out = 40; % Initial guess for outlet temperature (°C)
m_dot = 0.05;
Cp = 4183;
T_in = 24;
G = 912.286;
A_s = 1.85;
alpha_s = 0.96;

T_out_result = calculate_T_out(T_out, Ts, Tc, Ta, m_dot, Cp, T_in, G, A_s, alpha_s, V, L_c, air_properties, epsilon_c, sigma, error_boundary, L_s, beta_deg, g, epsilon_s);

%% Work out Ts from equation 1

function [Ts_new] = inner(Tc, Ts, T_out, air_properties, L_s, beta_deg, g, epsilon_s, epsilon_c, sigma, V, L_c, Ta, m_dot, Cp, T_in, h_i, A_i, error_boundary)
    abs_error = inf; % Initialize error to a large value
    Ts_new = Ts;     % Start with initial guess

    while abs_error > error_boundary
        % Step 1: Compute Convective Heat Transfer Coefficients
        [~, ~, ~, ~, h_s_c] = heatTransferCalc(Ts_new, Tc, L_s, beta_deg, g, air_properties);
    
        % Step 2: Compute Radiative Heat Transfer Coefficient
        hr_s_c = radiativeHeatTransfer(Ts_new, Tc, epsilon_s, epsilon_c, sigma);
     
        % Step 3: Compute Convective Heat Transfer Coefficient (Cover to Ambient)
        [~, ~, h_c_a] = calculateHeatTransfer(V, L_c, air_properties, Tc, Ta);
    
        % Step 4: Compute Radiative Heat Transfer Coefficient (Cover to Ambient)
        hr_c_a = calculateRadiativeHeatTransfer(Ta, Tc, epsilon_c, sigma);
    
        % Step 5: Compute Overall Heat Transfer Coefficient
        U_new = ((1 / (h_s_c + hr_s_c)) + (1 / (h_c_a + hr_c_a)))^(-1);
    
        % Step 6: Calculate new Ts using Equation 1
        Ts_old = Ts_new;
        % Compute log mean temperature difference (Delta T_LM)
        if abs(Ts_old - T_out) < 1e-6 || abs(Ts_old - T_in) < 1e-6
            % Avoid division by zero in log
            Delta_T_LM = (Ts_old - T_in + Ts_old - T_out) / 2; % Approximate as arithmetic mean
        else
            Delta_T_LM = ((Ts_old - T_in) - (Ts_old - T_out)) / log((Ts_old - T_in) / (Ts_old - T_out));
        end
        
        Ts_new = (m_dot * Cp * (T_out - T_in) / (h_i * A_i) + T_in + T_out) / 2; % Simplified initial update
        % More accurate: Solve iteratively by adjusting Ts
        for iter = 1:10 % Inner iteration to refine Ts
            if abs(Ts_new - T_out) < 1e-6 || abs(Ts_new - T_in) < 1e-6
                Delta_T_LM = (Ts_new - T_in + Ts_new - T_out) / 2;
            else
                Delta_T_LM = ((Ts_new - T_in) - (Ts_new - T_out)) / log((Ts_new - T_in) / (Ts_new - T_out));
            end
            Ts_new = (m_dot * Cp * (T_out - T_in) / (h_i * A_i) + T_in + T_out) / 2 + Delta_T_LM;
        end
    
        % Step 7: Compute absolute error
        abs_error = abs(Ts_new - Ts_old);
    
        % Step 8: Display progress
        fprintf('Current Ts = %.6f °C, Error = %.6f\n', Ts_new, abs_error);
    end
end

h_i = 1160.046918;
A_i = 7 * pi * 16.4e-3 * 1.8;
Ts_result = inner(Tc, Ts, T_out_result, air_properties, L_s, beta_deg, g, epsilon_s, epsilon_c, sigma, V, L_c, Ta, m_dot, Cp, T_in, h_i, A_i, error_boundary);

%% Work out Tc from equation 3

function [Tc_new] = outer(Tc, Ts, air_properties, L_s, beta_deg, g, epsilon_s, epsilon_c, sigma, V, L_c, Ta, A_s, error_boundary)
    abs_error = inf; % Initialize error to a large value
    Tc_new = Tc;     % Start with initial guess

    while abs_error > error_boundary
        % Step 1: Compute Convective Heat Transfer Coefficients
        [~, ~, ~, ~, h_s_c] = heatTransferCalc(Tc_new, Ts, L_s, beta_deg, g, air_properties);
    
        % Step 2: Compute Radiative Heat Transfer Coefficient
        hr_s_c = radiativeHeatTransfer(Tc_new, Ts, epsilon_s, epsilon_c, sigma);
     
        % Step 3: Compute Convective Heat Transfer Coefficient (Cover to Ambient)
        [~, ~, h_c_a] = calculateHeatTransfer(V, L_c, air_properties, Ts, Ta);
    
        % Step 4: Compute Radiative Heat Transfer Coefficient (Cover to Ambient)
        hr_c_a = calculateRadiativeHeatTransfer(Ta, Ts, epsilon_c, sigma);
    
        % Step 5: Compute Overall Heat Transfer Coefficient
        U_new = ((1 / (h_s_c + hr_s_c)) + (1 / (h_c_a + hr_c_a)))^(-1);
    
        % Step 6: Calculate new Tc using Equation 3
        Tc_old = Tc_new;
        % Equation 3: h_s_c * A_s * (Ts - Tc) = U * A_s * (Tc - Ta)
        % Rearrange: h_s_c * (Ts - Tc) = U * (Tc - Ta)
        % h_s_c * Ts - h_s_c * Tc = U * Tc - U * Ta
        % h_s_c * Ts + U * Ta = Tc * (h_s_c + U)
        % Tc = (h_s_c * Ts + U * Ta) / (h_s_c + U)
        Tc_new = ((h_s_c + hr_s_c) * Ts + U_new * Ta) / ((h_s_c + hr_s_c) + U_new);
    
        % Step 7: Compute absolute error
        abs_error = abs(Tc_new - Tc_old);
    
        % Step 8: Display progress
        fprintf('Current Tc = %.6f °C, Error = %.6f\n', Tc_new, abs_error);
    end
end

Tc_result = outer(Tc, Ts_result, air_properties, L_s, beta_deg, g, epsilon_s, epsilon_c, sigma, V, L_c, Ta, A_s, error_boundary);

%% Question 2: Convective Heat Transfer (Surface to Tube)
function [Tm, Gr, Ra, Nu_2, hc_s_c] = heatTransferCalc(Ts, Tc, L_s, beta_deg, g, air_properties)
   
    % Calculate mean temperature
    Tm = (Ts + Tc) / 2;    % Mean temperature in °C
    Tm_K = Tm + 273;       % Convert to Kelvin
    [k_m, nu_m, Pr_m] = interpolation(Tm, air_properties);
    % Grashof Number calculation
    beta = 1 / Tm_K; % Volumetric thermal expansion coefficient (1/K)
    Gr = (g * beta * (Ts - Tc) * L_s^3) / (nu_m^2);
    % Convert inclination angle to radians
    beta_rad = deg2rad(beta_deg);
    % Rayleigh Number
    Ra = Gr * Pr_m;
    % Corrected Holland Equation for Nusselt number
    Nu_2 = 1 + 1.44 * (1 - (1708 * (sin(1.8 * beta_rad)^1.6)) / (Ra * cos(beta_rad))) * ...
         max(0, 1 - (1708 / (Ra * cos(beta_rad)))) + ...
         max(0, (Ra * cos(beta_rad) / 5830)^(1/3) - 1);
    % Average heat transfer coefficient
    hc_s_c = (Nu_2 * k_m) / L_s;
    % Display results
    fprintf('Convection Heat Transfer Coefficient (hc_s_c): %.6f W/(m^2·K)\n', hc_s_c);
end
 
%% Question 3: Radiative Heat Transfer (Surface to Tube)
function hr_s_c = radiativeHeatTransfer(Ts, Tc, epsilon_s, epsilon_c, sigma)

    Ts_K = Ts + 273; % Surface temperature in Kelvin
    Tc_K = Tc + 273; % Ambient temperature in Kelvin
    % Calculate radiative heat transfer coefficient (hr_s_c)
    hr_s_c = (sigma) * ((Ts_K)^2 + (Tc_K)^2) * (Ts_K + Tc_K) / ((1/epsilon_s) + (1/epsilon_c) - 1);
    % Display result
    fprintf('Radiative Heat Transfer Coefficient (hr_s_c): %.6f W/(m^2·K)\n', hr_s_c);
end
 
%% Question 4: Convective Heat Transfer (Tube to Ambient)
function [Re, Nu_4, hc_c_a] = calculateHeatTransfer(V, L_c, air_properties, Tc, Ta)
   
    Tm = (Tc + Ta) / 2;
    [k_a, nu_a, Pr_a] = interpolation(Tm, air_properties);
    % Calculate Reynolds number (Re)
    Re = (V * L_c) / nu_a;
    % Calculate Nusselt number (Nu) using the correlation Nu = 0.664 * Re^(1/2) * Pr^(1/3)
    Nu_4 = 0.664 * (Re^(1/2)) * (Pr_a^(1/3));
    % Calculate convective heat transfer coefficient (h_c-a)
    hc_c_a = (Nu_4 * k_a) / L_c;
    % Display the results
    fprintf('Convective Heat Transfer Coefficient (h_c-a) = %.6f W/m^2·K\n', hc_c_a);
end
 
%% Question 5: Radiative Heat Transfer (Tube to Ambient)
function hr_c_a = calculateRadiativeHeatTransfer(Ta, Tc, epsilon_c, sigma)
    
    Ta_5 = Ta + 273;         % Ambient temperature in Kelvin
    Tc_5 = Tc + 273;         % Tube surface temperature in Kelvin
    % Calculate radiative heat transfer coefficient (h_r,c-a)
    hr_c_a = epsilon_c * sigma * (Ta_5^2 + Tc_5^2) * (Ta_5 + Tc_5);
    % Display the result
    fprintf('Radiative Heat Transfer Coefficient (h_r,c-a) = %.6f W/m^2·K\n', hr_c_a);
end
 
%% Interpolation Function

function [T1, T2] = find_temp_bounds(Tm, air_properties)
    temperatures = air_properties(:,1);
    temperatures = sort(temperatures);
 
    if Tm < temperatures(1) || Tm > temperatures(end)
        error('Tm = %.6f is out of the table range (%.6f to %.6f).', Tm, temperatures(1), temperatures(end));
    end
 
    for i = 1:length(temperatures) - 1
        if temperatures(i) <= Tm && Tm <= temperatures(i + 1)
            T1 = temperatures(i);
            T2 = temperatures(i + 1);
            return;
        end
    end
 
    error('Unexpected error in find_temperature_bounds. Check input values.');
end
 
function [k1, k2, nu1, nu2, Pr1, Pr2] = air_properties_calc(T1, T2, air_properties)
    idx1 = find(air_properties(:,1) == T1, 1);
    idx2 = find(air_properties(:,1) == T2, 1);
 
    if isempty(idx1) || isempty(idx2)
        error('Temperature values not found in air properties table.');
    end
 
    k1 = air_properties(idx1, 2);
    k2 = air_properties(idx2, 2);
    nu1 = air_properties(idx1, 3);
    nu2 = air_properties(idx2, 3);
    Pr1 = air_properties(idx1, 4);
    Pr2 = air_properties(idx2, 4);
end
 
function [k, nu, Pr] = interpolation(Tm, air_properties)
    temperatures = air_properties(:,1);
    Tmin = min(temperatures);
    Tmax = max(temperatures);

    Tm = max(min(Tm, Tmax), Tmin);

    [T1, T2] = find_temp_bounds(Tm, air_properties);
    [k1, k2, nu1, nu2, Pr1, Pr2] = air_properties_calc(T1, T2, air_properties);
 
    k = k1 + (k2 - k1) * ((Tm - T1) / (T2 - T1));
    nu = nu1 + (nu2 - nu1) * ((Tm - T1) / (T2 - T1));
    Pr = Pr1 + (Pr2 - Pr1) * ((Tm - T1) / (T2 - T1));
end

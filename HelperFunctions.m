%% Question 3 function
function t_220 = calculate_time_to_220_degrees()
    % Function to calculate the time for the gas temperature to reach 220°C
    % using the lumped capacitance method with a linearly increasing surrounding temperature
    
    % Given parameters
    r = (0.25e-3) / (2 * 3);        % Radius of the spherical sensor (m)
    h = 60;                         % Convective heat transfer coefficient (W/m^2·K)
    k = 70;                         % Thermal conductivity of the gas (W/m·K)
    rho_cp = 21452;                 % Volumetric heat capacity (J/m^3·K)
    beta = 200;                     % Rate of temperature increase (°C/s)
    T_ini = 20;                     % Initial surrounding temperature (°C)
    T_target = 220;                 % Target temperature (°C)
    
    % Step 1: Calculate characteristic length (l_c = r for this problem)
    l_c = r;  % Characteristic length (m)
    fprintf('Characteristic Length (l_c) = %.2e m\n', l_c);
    
    % Step 2: Calculate Biot number (Bi)
    Bi = (h * l_c) / k;
    fprintf('Biot Number (Bi) = %.2e\n', Bi);
    
    % Check if lumped capacitance is valid (Bi < 0.1)
    if Bi >= 0.1
        error('Biot number (Bi = %.2e) is not less than 0.1. Lumped capacitance method is not valid.', Bi);
    end
    
    % Step 3: Calculate lumped capacitance time constant (tau_lumped)
    V = (4/3) * pi * r^3;          % Volume of the sphere (m^3)
    A = 4 * pi * r^2;              % Surface area of the sphere (m^2)
    tau_lumped = (rho_cp * V) / (h * A);
    fprintf('Lumped Capacitance Time Constant (tau_lumped) = %.3f s\n', tau_lumped);
    
    % Step 4: Define the temperature equation and solve for t
    % T(t) = beta * tau_lumped * (exp(-t/tau_lumped) - 1) + T_ini + beta * t
    f = @(t) beta * tau_lumped * (exp(-t/tau_lumped) - 1) + T_ini + beta * t - T_target;
    
    % Use fzero to solve for t when T = 220°C
    t_initial_guess = 1;  % Initial guess for t (s)
    t_220 = fzero(f, t_initial_guess);
    
    % Display the result
    fprintf('Time to reach 220°C: t_220 = %.3f s\n', t_220);
    fprintf('Rounded to 1 decimal place: t_220 = %.1f s\n', round(t_220, 1));
end

%% Solve for the time when the gas temperature reaches 320°C

% Define the equation: t^2 - 2.96*t - 44.14 = 5.56*exp(-t)
f = @(t) t^2 - 2.96*t - 44.14 - 5.56*exp(-t);

% Use fzero to solve for t
t_initial_guess = 8;  % Initial guess for t (s), close to the expected solution
t_320 = fzero(f, t_initial_guess);

% Display the result
fprintf('Time to reach 320°C: t_320 = %.4f s\n', t_320);

%% Interpolating

% Define the known data points from Table A-9
temp = [35 40]; % Temperatures in °C

% Properties at 30°C and 40°C
density = [1.145 1.127];          % kg/m³
heat_capacity = [1007 1007];      % J/kg·K
conductivity = [0.02625 0.02662]; % W/m·K
diffusivity = [2.277e-5 2.346e-5]; % m²/s
dynamic_viscosity = [1.895e-5 1.918e-5]; % kg/m·s
kinematic_viscosity = [1.655e-5 1.702e-5]; % m²/s
prandtl_number = [0.7268 0.7255]; % Pr

% Target temperaturem
T_target = 35.5; % °C

% Perform linear interpolation for each property
rho_interp = interp1(temp, density, T_target, 'linear');
cp_interp = interp1(temp, heat_capacity, T_target, 'linear');
k_interp = interp1(temp, conductivity, T_target, 'linear');
alpha_interp = interp1(temp, diffusivity, T_target, 'linear');
mu_interp = interp1(temp, dynamic_viscosity, T_target, 'linear');
nu_interp = interp1(temp, kinematic_viscosity, T_target, 'linear');
Pr_interp = interp1(temp, prandtl_number, T_target, 'linear');

% Display results
fprintf('Interpolated Air Properties at %.1f°C and 1 atm pressure:\n', T_target);
fprintf('Density (ρ): %.3f kg/m³\n', rho_interp);
fprintf('Specific Heat Capacity (c_p): %.0f J/kg·K\n', cp_interp);
fprintf('Thermal Conductivity (k): %.5f W/m·K\n', k_interp);
fprintf('Thermal Diffusivity (α): %.5e m²/s\n', alpha_interp);
fprintf('Dynamic Viscosity (μ): %.5e kg/m·s\n', mu_interp);
fprintf('Kinematic Viscosity (ν): %.5e m²/s\n', nu_interp);
fprintf('Prandtl Number (Pr): %.4f\n', Pr_interp);

%% Interpolating

% Define the known data points from Table A-9
temp = [20 25]; % Temperatures in °C

% Properties at 20°C and 25°C
density = [1.204 1.184];          % kg/m³
heat_capacity = [1007 1007];      % J/kg·K
conductivity = [0.02514 0.02551]; % W/m·K
diffusivity = [2.074e-5 2.141e-5]; % m²/s
dynamic_viscosity = [1.825e-5 1.849e-5]; % kg/m·s
kinematic_viscosity = [1.516e-5 1.562e-5]; % m²/s
prandtl_number = [0.7309 0.7296]; % Pr

% Target temperature
T_target = 21; % °C

% Perform linear interpolation for each property
rho_interp = interp1(temp, density, T_target, 'linear');
cp_interp = interp1(temp, heat_capacity, T_target, 'linear');
k_interp = interp1(temp, conductivity, T_target, 'linear');
alpha_interp = interp1(temp, diffusivity, T_target, 'linear');
mu_interp = interp1(temp, dynamic_viscosity, T_target, 'linear');
nu_interp = interp1(temp, kinematic_viscosity, T_target, 'linear');
Pr_interp = interp1(temp, prandtl_number, T_target, 'linear');

% Display results
fprintf('Interpolated Air Properties at %.1f°C and 1 atm pressure:\n', T_target);
fprintf('Density (ρ): %.3f kg/m³\n', rho_interp);
fprintf('Specific Heat Capacity (c_p): %.0f J/kg·K\n', cp_interp);
fprintf('Thermal Conductivity (k): %.5f W/m·K\n', k_interp);
fprintf('Thermal Diffusivity (α): %.5e m²/s\n', alpha_interp);
fprintf('Dynamic Viscosity (μ): %.5e kg/m·s\n', mu_interp);
fprintf('Kinematic Viscosity (ν): %.5e m²/s\n', nu_interp);
fprintf('Prandtl Number (Pr): %.4f\n', Pr_interp);

%% Calculating Nusslet


function Nu = calculateNu(Ra, beta)
    % Function to calculate Nusselt number (Nu) based on Rayleigh number (Ra) and angle (beta)
    
    % Convert beta to radians if it's in degrees
    beta_rad = deg2rad(beta);
    
    % Calculate Nu in a single expression
    Nu = 1 + 1.44 * (1 - (1708 * (sin(1.8 * beta_rad)^1.6)) / (Ra * cos(beta_rad))) * ...
         max(0, 1 - (1708 / (Ra * cos(beta_rad)))) + ...
         max(0, (Ra * cos(beta_rad) / 5830)^(1/3) - 1);
end

% Input values
beta = 10;      % Angle in degrees
Ra = 1.043e5;   % Rayleigh number
Nu = calculateNu(Ra, beta);
disp(['Nusselt number (Nu) = ', num2str(Nu)]);

%% Interpolation

% Data from the table for 20°C and 25°C
T = [20, 25]; % Temperatures in °C

% Values at 20°C and 25°C
rho = [998.2, 997.1];              % Density (kg/m^3)
mu = [0.001002, 0.00089005];       % Dynamic viscosity (Pa·s)
Cp = [4183, 4183];                 % Specific heat (J/kg·K)
k = [0.5861, 0.5948];              % Thermal conductivity (W/m·K)
Pr = [7.152, 6.187];               % Prandtl number
beta = [0.000209, 0.000294];       % Thermal expansion coefficient (1/K)
c = [1483, 1501];                  % Sound speed (m/s)
sigma = [0.07273, 0.07119];        % Surface tension (N/m)

% Target temperature
T_interp = 24;

% Linear interpolation
rho_24 = interp1(T, rho, T_interp);
mu_24 = interp1(T, mu, T_interp);
Cp_24 = interp1(T, Cp, T_interp);
k_24 = interp1(T, k, T_interp);
Pr_24 = interp1(T, Pr, T_interp);
beta_24 = interp1(T, beta, T_interp);
c_24 = interp1(T, c, T_interp);
sigma_24 = interp1(T, sigma, T_interp);

% Display interpolated results
fprintf('Interpolated thermophysical properties of water at 24°C:\n');
fprintf('Density (ρ): %.4f kg/m³\n', rho_24);
fprintf('Dynamic viscosity (μ): %.8f Pa·s\n', mu_24);
fprintf('Specific heat (Cp): %.2f J/kg·K\n', Cp_24);
fprintf('Thermal conductivity (k): %.4f W/m·K\n', k_24);
fprintf('Prandtl number (Pr): %.4f\n', Pr_24);
fprintf('Thermal expansion coefficient (β): %.9f 1/K\n', beta_24);
fprintf('Sound speed (c): %.2f m/s\n', c_24);
fprintf('Surface tension (σ): %.5f N/m\n', sigma_24);

%% Interpolation

% Given temperatures and corresponding properties
T1 = 25; % °C
T2 = 30; % °C
T = 27.183; % °C (temperature to interpolate at)

% Properties at T1 (25°C) and T2 (30°C)
% Density (kg/m³)
rho1 = 997.1;
rho2 = 995.7;

% Dynamic Viscosity (Pa·s)
mu1 = 0.0008905;
mu2 = 0.0007977;

% Specific Heat (J/kg·K)
cp1 = 4183;
cp2 = 4183;

% Thermal Conductivity (W/m·K)
k1 = 0.5948;
k2 = 0.603;

% Prandtl Number
Pr1 = 6.263;
Pr2 = 5.534;

% Volume Expansion Coefficient (K⁻¹)
beta1 = 0.0002594;
beta2 = 0.0003051;

% Speed of Sound (m/s)
c1 = 1497;
c2 = 1509;

% Surface Tension (N/m)
sigma1 = 0.07197;
sigma2 = 0.07119;

% Linear interpolation formula: y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
% Interpolating each property at T = 27.183°C
rho = rho1 + (rho2 - rho1) * (T - T1) / (T2 - T1);
mu = mu1 + (mu2 - mu1) * (T - T1) / (T2 - T1);
cp = cp1 + (cp2 - cp1) * (T - T1) / (T2 - T1);
k = k1 + (k2 - k1) * (T - T1) / (T2 - T1);
Pr = Pr1 + (Pr2 - Pr1) * (T - T1) / (T2 - T1);
beta = beta1 + (beta2 - beta1) * (T - T1) / (T2 - T1);
c = c1 + (c2 - c1) * (T - T1) / (T2 - T1);
sigma = sigma1 + (sigma2 - sigma1) * (T - T1) / (T2 - T1);

% Display the results
fprintf('Interpolated thermophysical properties of water at T = %.3f°C:\n', T);
fprintf('Density (ρ): %.2f kg/m³\n', rho);
fprintf('Dynamic Viscosity (μ): %.7f Pa·s\n', mu);
fprintf('Specific Heat (cₚ): %.1f J/kg·K\n', cp);
fprintf('Thermal Conductivity (k): %.4f W/m·K\n', k);
fprintf('Prandtl Number (Pr): %.3f\n', Pr);
fprintf('Volume Expansion Coefficient (β): %.7f K⁻¹\n', beta);
fprintf('Speed of Sound (c): %.1f m/s\n', c);
fprintf('Surface Tension (σ): %.5f N/m\n', sigma);

%% Working out h_i

% Calculate the convective heat transfer coefficient inside the tube (h_i)

% Given parameters
D = 16.4e-3;        % Diameter of the tube (m)
m_dot = 0.05;       % Mass flow rate (kg/s)
T_water = 24;       % Water temperature (°C), assumed as T_in

% Water properties at 24°C (interpolated between 20°C and 25°C from a typical table)
% Using approximate values for water (since the code uses air properties, we define water properties here)
mu = 9.1244e-4;     % Dynamic viscosity of water at 24°C (kg/m·s)
Pr = 6.3800;        % Prandtl number of water at 24°C
k = 0.5951;         % Thermal conductivity of water at 24°C (W/m·K)

% Calculate velocity (V)
A = (pi * D^2) / 4; % Cross-sectional area of the tube (m^2)
V = m_dot / (A * 997.32); % Velocity (m/s), using water density at 24°C (997.32 kg/m^3)

% Calculate Reynolds number (Re)
Re = (997.32 * V * D) / mu;
fprintf('Reynolds Number (Re) = %.2f\n', Re);

% Calculate Nusselt number (Nu) using Dittus-Boelter correlation: Nu = 0.023 * Re^0.8 * Pr^0.3
Nu = 0.023 * Re^0.8 * Pr^0.3;
fprintf('Nusselt Number (Nu) = %.2f\n', Nu);

% Calculate convective heat transfer coefficient (h_i)
h_i = (Nu * k) / D;
fprintf('Convective Heat Transfer Coefficient (h_i) = %.6f W/m^2·K\n', h_i);

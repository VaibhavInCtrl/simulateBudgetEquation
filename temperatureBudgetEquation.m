function result = temperatureBudgetEquation()
    % Gradient calculation function
    function gradient = get_gradient(x, y)
        p = polyfit(x, y, 1);
        gradient = p(1);
    end

    % Function to calculate penetrative shortwave radiation (Qpen)
    h = 50; % Mixed layer depth (m)
    QSWR = 300; % Shortwave radiation flux (W/m^2)
    R = 0.58; % Reflection coefficient for shortwave radiation
    zeta1 = 0.35; % Attenuation coefficient for the first layer (1/m)
    zeta2 = 23; % Attenuation coefficient for the second layer (1/m)
    function penetration = Qpen(QSWR, h, R, zeta1, zeta2)
        penetration = QSWR * (R * exp(-h/zeta1) + (1 - R) * exp(-h/zeta2));
    end

    % Parameters for temperature advection calculation
    tau_x = 0.1; % Wind stress in x-direction (N/m^2)
    tau_y = 0.1; % Wind stress in y-direction (N/m^2)
    f = 1e-4; % Coriolis parameter (1/s), typical for mid-latitudes
    rho = 1025; % Density of seawater (kg/m^3)
    additional_currents = [0.05, 0.05]; % Other currents (u_other, v_other) in m/s
    function horizontal_advection = calculate_temperature_advection(dTdx, dTdy)
        u_e = tau_y / (rho * f * h);
        v_e = -tau_x / (rho * f * h);
        u_total = u_e + additional_currents(1);
        v_total = v_e + additional_currents(2);
        horizontal_advection = u_total * dTdx + v_total * dTdy;
    end

    % Main function to compute temperature residual
    function residual = temperature_residual(T, dTdx, dTdy, wh, Th, QSWR, dTdt)
        penetration = Qpen(QSWR, h, R, zeta1, zeta2);
        H = double(wh > 0); % Heaviside function equivalent
        Qnet = 100; % Net surface heat flux (W/m^2)
        cp = 3985; % Specific heat capacity at constant pressure (J/(kg*K))
        Kz = 1e-5; % Vertical diffusivity (m^2/s)
        heating = (Qnet - penetration) / (rho * cp * h);
        advection = -calculate_temperature_advection(dTdx, dTdy);
        entrainment = -(H * (wh + dhdt) * ((T - Th) / h));
        diffusion = -(Kz / h) * dTdz;
        residual = dTdt - heating - advection - entrainment - diffusion;
    end

    % Example data and calculations
    h_array = [50,55,60];
    t_array = [10,13,17];
    dhdt = get_gradient(h_array, t_array);
    T_array = [35,37,40];
    z_array = [2,3,5];
    dTdz = get_gradient(T_array, z_array);
    T = 20;
    Th = 19;
    x_array = [3,5,8];
    dTdx = get_gradient(T_array, x_array);
    y_array = [5,6,7];
    dTdy = get_gradient(T_array, y_array);
    wh = 0.01;
    dTdt = get_gradient(T_array, t_array);
    residual = temperature_residual(T, dTdx, dTdy, wh, Th, QSWR, dTdt);
    fprintf('Temperature Budget Residual: %.6f K/s\n', residual);
end

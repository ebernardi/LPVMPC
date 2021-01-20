function mu = membership(x, rho1_min, rho1_max, rho2_min, rho2_max)
    % Estados
    % x(1) => Tf
    % x(2) => Tp

    %% Parameters
    mum = 1100;       % [kg/m³]
    muf = 1000;       % [kg/m³]
    Ai = 0.013;       % [m²]
    Ae = 0.0038;      % [m²]
    di = 0.04;        % [m]
    de = 0.07;        % [m]
    h0 = 11;
    hi = 800;
    Cm = 440;         % [J/(kgC)]
    Cf = 4018;        % [J/(kgC)]
    nhu = 0.43*8.5;
    
    rho1 = di*pi*hi*(1-exp(-x(1)/600))/(1-exp(-1));
    rho2 = (1-exp(-x(2)/300))/((1-exp(-1))*Ai);
    
    % Membresía rho1
    M1_1 = (rho1_max-rho1)/(rho1_max-rho1_min);
    M1_2 = (rho1-rho1_min)/(rho1_max-rho1_min);
    
    % Membresía rho2
    M2_1 = (rho2_max-rho2)/(rho2_max-rho2_min);
    M2_2 = (rho2-rho2_min)/(rho2_max-rho2_min);        

    mu1 = M1_1*M2_1;
    mu2 = M1_1*M2_2;
    mu3 = M1_2*M2_1;
    mu4 = M1_2*M2_2;
    
    mu = [mu1 mu2 mu3 mu4];
end

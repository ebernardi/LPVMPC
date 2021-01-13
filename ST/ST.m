%% Solar Collector non-linear model
function out = ST(t, x, u, w)
	%
	% Solar Collector function
	%
	% x: State vector [Tp Tf]
	% u: Input vector [u]
    % w: Disturbance vector [Irr Te]
	% out: Output vector [dTp dTf]

	% States
	Tp = x(1);
	Tf = x(2);
	
    % Input
    u = u(1);
    
    % Disturbances
    Irr = w(1);
    Te = w(2);
		
	% Parameters
    mum = 1100;       % [kg/m³]
    muf = 1000;         % [kg/m³]
    Ai = 0.013;           % [m²]
    Ae = 0.0038;        % [m²]
    di = 0.04;             % [m]
    de = 0.07;            % [m]
    h0 = 11;
    hi = 800;
    Cm = 440;           % [J/(kgC)]
    Cf = 4018;           % [J/(kgC)]
    nhu = 0.43*8.5;

    rho_1 = di*pi*hi*(1-exp(-Tp/600))/(1-exp(-1));
    rho_2 = (1-exp(-Tf/300))/((1-exp(-1))*Ai);
    
	% Energy balance
	dTp = (de*pi*nhu*Irr - de*pi*h0*(Tp - Te) - rho_1*(Tp - Tf))/(mum*Cm*Ae);
	dTf = -u*rho_2 + rho_1*(Tp - Tf)/(muf*Cf*Ai);
	
 	out = [dTp; dTf];
end
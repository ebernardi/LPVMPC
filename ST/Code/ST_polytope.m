%% ST Model
syms Uf TF TP I TE Rho1 Rho2
sys = struct;

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

%% Non-Linear model
rho1 = di*pi*hi*(1-exp(-TP/600))/(1-exp(-1));
rho2 = (1-exp(-TF/300))/((1-exp(-1))*Ai);
system = [(de*pi*nhu*I-de*pi*h0*(TP-TE)-Rho1*(TP-TF))/(mum*Cm*Ae);
                  -Uf*Rho2+(Rho1*(TP-TF)/(muf*Cf*Ai))];

states = [TP TF];
outputs = [TP TF];
disturbances = [I TE];
parameters = [Rho1 Rho2];
inputs = Uf;
nx = length(states); nu = length(inputs);  ny = 1; nd = length(disturbances);
C = [0 1];       % Output matrix
D = zeros(ny, nu);    % Input/Output matrix

%% Linealization
% Symbolic matrices
A_sym = jacobian(system, states);
B_sym = jacobian(system, inputs);
Bw_sym = jacobian(system, disturbances);

for i = 1:M
	if M == 4
        switch i
            case 1
                % System 1 (TP = Tp_min [K]; TF = Tf_min [K])
                Tp = Tp_min;		   % [K] Plate temperature
                Tf = Tf_min;             % [K] Fluid temperature
            case 2
                % System 2 (TP = Tp_min [K]; TF = Tf_max [K])
                Tp = Tp_min;			% [K] Plate temperature
                Tf = Tf_max;             % [K] Fluid temperature
            case 3
                % System 3 (TP = Tp_max [K]; TF = Tf_min [K])
                Tp = Tp_max;		  % [K] Plate temperature
                Tf = Tf_min;             % [K] Fluid temperature
            case 4
                % System 4 (TP = Tp_max [K]; TF = Tf_max [K])
                Tp = Tp_max;	 	   % [K] Plate temperature
                Tf = Tf_max;             % [K] Fluid temperature
            otherwise
                break;
        end
    elseif M == 9
        switch i
            case 1
                % System 1 (TP = Tp_min [K]; TF = Tf_min [K])
                Tp = Tp_min;		   % [K] Plate temperature
                Tf = Tf_min;             % [K] Fluid temperature
            case 2
                % System 2 (TP = Tp_min [K]; TF = Tf_mid [K])
                Tp = Tp_min;			% [K] Plate temperature
                Tf = Tf_mid;             % [K] Fluid temperature
            case 3
                % System 3 (TP = Tp_min [K]; TF = Tf_max [K])
                Tp = Tp_min;		  % [K] Plate temperature
                Tf = Tf_max;             % [K] Fluid temperature
            case 4
                % System 4 (TP = Tp_mid [K]; TF = Tf_min [K])
                Tp = Tp_mid;	 	   % [K] Plate temperature
                Tf = Tf_min;             % [K] Fluid temperature
            case 5
                % System 5 (TP = Tp_mid [K]; TF = Tf_mid [K])
                Tp = Tp_mid;		   % [K] Plate temperature
                Tf = Tf_mid;             % [K] Fluid temperature
            case 6
                % System 6 (TP = Tp_mid [K]; TF = Tf_max [K])
                Tp = Tp_mid;			% [K] Plate temperature
                Tf = Tf_max;             % [K] Fluid temperature
            case 7
                % System 7 (TP = Tp_max [K]; TF = Tf_min [K])
                Tp = Tp_max;		  % [K] Plate temperature
                Tf = Tf_min;             % [K] Fluid temperature
            case 8
                % System 8 (TP = Tp_max [K]; TF = Tf_mid [K])
                Tp = Tp_max;	 	   % [K] Plate temperature
                Tf = Tf_mid;             % [K] Fluid temperature
            case 9
                % System 9 (TP = Tp_max [K]; TF = Tf_max [K])
                Tp = Tp_max;	 	   % [K] Plate temperature
                Tf = Tf_max;             % [K] Fluid temperature
            otherwise
                break;
        end
	end
    % Scheduling parameters
    rho1 = di*pi*hi*(1-exp(-Tp/600))/(1-exp(-1));
    rho2 = (1-exp(-Tf/300))/((1-exp(-1))*Ai);
    
    rho1_min = di*pi*hi*(1-exp(-Tp_min/600))/(1-exp(-1));
    rho1_max = di*pi*hi*(1-exp(-Tp_max/600))/(1-exp(-1));
    rho2_min = (1-exp(-Tf_min/300))/((1-exp(-1))*Ai);
    rho2_max = (1-exp(-Tf_max/300))/((1-exp(-1))*Ai);
    
    % Disturbances
    Winit = [Irr; Te];
    
    % States
    Xinit = [Tp; Tf];

    % Fluid flow rate
    u = solve(-Uf*rho2 + rho1*(Tp - Tf)/(muf*Cf*Ai)==0);
    u = double(u);
    
    % Inputs
    Uinit = u;

    % Linear systems matrices
    A = subs(A_sym, {Rho1, Rho2}, {rho1, rho2});
    B = subs(B_sym, {Rho1, Rho2}, {rho1, rho2});
    Bw = subs(Bw_sym, {I, TE}, {Irr, Te});
    A = double(A);
    B = double(B);
    Bw = double(Bw);
    
    f = subs(system, {TP, TF, Rho1, Rho2, Uf, I, TE}, {Tp, Tf, rho1, rho2, u, Irr, Te});
    f = double(f);

    % Constant term
    delta = f - (A*Xinit+B*Uinit+Bw*Winit);

	% Continuous system
    sys(i).A = A; sys(i).B = B; sys(i).delta = delta;
    sys(i).C = C; sys(i).D = D;
    
    % Steady state
    sys(i).Xinit = Xinit; sys(i).Uinit = Uinit; sys(i).Winit = Winit;
    
    % Euler discretization method
    sys(i).Ad = (A*Ts) + eye(nx); sys(i).Bd = B*Ts; sys(i).Bwd = Bw*Ts; sys(i).deltad = delta*Ts;
    
    % Common matrices
    sys(i).Cd = C; sys(i).Dd = D;
end
%% LPV MPC on ST
clc; clear all; yalmip('clear');
close all;

load('Datos Acurex 2007-06-18.mat')

% ODE options 'RelTol', 1e-6, 'AbsTol', 1e-6
options = odeset ('RelTol', 1e-6, 'AbsTol', 1e-6, ...
	'NormControl', 'on', 'InitialStep', 1.0e-2, 'MaxStep', 1.0);

% Ts = 3;                     % Sample time [seg]
% t = 0;                       % Start time
% Time = 3.6e3;          % Simulation end time 
% Tsim = 0:Ts:Time;
% Nsim = length(Tsim);
% 
% % Disturbances
% Irr = IC151(1, 2);
% Te = TA075(1, 2);
% w0 = [Irr; Te];              % Init disturbance [Irr; Te]
% 
% %% Polytope model and sets
% % This section is commented to reduce simulation time (using pre-calculated matrices)
% Tp_max = 600;
% Tp_min = 0.0377;% @ rho1 = 0.01
% Tf_max = 300;
% Tf_min = 0.0247; % @ rho2 = 0.01
% 
% N = 2;                       % Number of parameters
% L = 2;                        % Linearization points per parameter
% M = L^N;                  % Number of models
% 
% run ST_polytope
% 
% %% Controller parameters
% N_MPC = 10;
% N_MHE = N_MPC;
% 
% % Box-type constraints
% xmin = [Tp_min; Tf_min];
% xmax = [Tp_max; Tf_max];
% umin = 0;
% umax = 0.35;
% Xset = Polyhedron('lb', xmin, 'ub', xmax);
% Uset = Polyhedron('lb', umin, 'ub', umax);
% 
% % Wheight matrices
% Qx = 1e4*sys(1).Cd'*sys(1).Cd;
% Ru = 1e1;
% Qe = 1e5;
% Qnux = 1;
% Qnumu =8e2;
% 
% % Set-point (Tf_sp = 97)
% Tp_sp = 109.92832;
% rho1 = di*pi*hi*(1-exp(-Tp_sp/600))/(1-exp(-1));
% syms Tf u
% Tf_sp = solve((-de*pi*h0/(mum*Cm*Ae)-rho1/(mum*Cm*Ae))*Tp_sp+rho1*Tf/(mum*Cm*Ae)+de*pi*nhu*Irr/(mum*Cm*Ae)+de*pi*h0*Te/(mum*Cm*Ae) == 0);
% Tf_sp = double(Tf_sp);
% Xsp = [Tp_sp; Tf_sp];
% rho2 = (1-exp(-Tf_sp/300))/((1-exp(-1))*Ai);
% Usp = solve(-u*rho2 + rho1*(Tp_sp - Tf_sp)/(muf*Cf*Ai)==0);
% Usp = double(Usp);
% 
% % |u| ≤ Usp
% u_abs = Usp; % umax-u_ss (can't handle asymmetric input constraints)
% 
% % Start point
% x0 = [300; 92.5];
% rho1 = di*pi*hi*(1-exp(-x0(1)/600))/(1-exp(-1));
% rho2 = (1-exp(-x0(2)/300))/((1-exp(-1))*Ai);
% syms u
% u0 = solve(-u*rho2 + rho1*(x0(1) - x0(2))/(muf*Cf*Ai)==0);
% u0 = double(u0);
% 
% % Terminal ingredients
% run terminalSetsCisneros
% % Save data
% save data.mat

load data

%% MHE
run MHE

%% MPC
run MPC

%% Simulation
% Vector initialization for plots
Obj = []; input = []; W = []; Y = C*x0; X = x0; Xd = x0; deltaU = zeros(1, Nsim);
X_MHE = repmat(x0, 1, N_MHE+1); U_MHE = repmat(u0, 1, N_MHE); 
mu_mhe = ones(M, Nsim)/M; mu_fuzzy = zeros(M, Nsim);
mu_MHE = mu_mhe(:, 1);

Predict = zeros(nd, N_MPC);
for k = 1:N_MPC
            Predict(:, k) = [IC151(k+1, 2); TA075(k+1, 2)];
end

disp('Iniciando...')
time_MHE = zeros(Nsim, 1) ;  % initilaize the elapsed times 
time_MPC = zeros(Nsim, 1) ;  % initilaize the elapsed times 
    for i = 1:Nsim
        tk = i*Ts;
        disp(i);
        
        % Update disturbance
        Irr = IC151(i, 2);
        Te = TA075(i, 2);
        w = [Irr Te]';
        
        t_tic = tic;   % to get time evaluated 
        
        [sol, diag] = mhe{X_MHE, U_MHE, Predict, mu_MHE};
        if diag
           msg = ['Infeasible MHE at t = ', num2str(i*Ts)];
           disp(msg)
           return;
        end
        mu_mhe(:, i) = sol;
        
        t_tic = toc(t_tic);         % get time elapsed
        time_MHE(i) = t_tic;   % store the time elapsed for the run
        mu_fuzzy(:, i) = membership(x0, rho1_min, rho1_max, rho2_min, rho2_max)';
        
        % Disturbance prediction
        for k = 1:N_MPC
            Predict(:, k) = [IC151(i+k+1, 2); TA075(i+k+1, 2)];
        end
        
        t_tic = tic;   % to get time evaluated 
        
        [sol, diag] = mpc{x0, Xsp, Usp, Predict, mu_mhe(:, i)};
        if diag
            msg = ['Infeasible MPC at t = ', num2str(i*Ts)];
            disp(msg)
            return;
        end
        umpc = sol{1}; 
        xpred = sol{2};
        obj = sol{3};
        
        t_tic = toc(t_tic) ;              % get time elapsed
        time_MPC(i) = t_tic;   % store the time elapsed for the run

        % Continuous-time simulation (reality)
        [tsim, x] = ode45(@ST, [0 Ts], x0, options, umpc, w);
        x0 = x(end, :)';                 % Next state

        X = [X x(2:end, :)'];          % State vector
        Xd = [Xd x(end, :)'];          % Discrete state vector
        Y = [Y sys(1).Cd*x(2:end, :)'];     % Output vector
        t = [t tk+tsim(2:end)'];    % Time vector
        input = [input; umpc];     % Input vector
        Obj = [Obj obj];                % Objective vector
        W = [W w];                        % Disturbance vector
        
        % Update MHE
        X_MHE = [X_MHE(:, 2:end) x0];
        U_MHE = [U_MHE(:, 2:end) umpc];
        mu_MHE = mu_mhe(:, i);
        
        deltaU(i) = abs(U_MHE(end)-U_MHE(end-1)); 

    end

%% Save and plot data
save run.mat
run enPlotST
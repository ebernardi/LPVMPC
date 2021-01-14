%% LPV MPC on ST
% clc; clear all; yalmip('clear');
close all;

load('Datos Acurex 2007-06-18.mat')

% ODE options 'RelTol', 1e-6, 'AbsTol', 1e-6
options = odeset ('RelTol', 1e-6, 'AbsTol', 1e-6, ...
	'NormControl', 'on', 'InitialStep', 1.0e-2, 'MaxStep', 1.0);

Ts = 3;                % Sample time [seg]
N = 8;                 % Prediction horizon
t = 0;                  % Start time
Time = 3.6e3;          % Simulation end time 
Tsim = 0:Ts:Time;
Nsim = length(Tsim);

% Disturbances
Irr = IC151(1, 2);
Te = TA075(1, 2);
w0 = [Irr; Te];              % Init disturbance [Irr; Te]

Tf_max = 300;
Tf_min = 0.0247; % @ rho2 = 0.01
Tp_max = 600;
Tp_min = 0.0377;% @ rho1 = 0.01

L = 2;  % 2 linearization points per parameter

run ST_polytope

% States
x0 = [300; 92.5];            % Start-point
% U
rho1 = di*pi*hi*(1-exp(-x0(1)/600))/(1-exp(-1));
rho2 = (1-exp(-x0(2)/300))/((1-exp(-1))*Ai);
syms u
u0 = solve(-u*rho2 + rho1*(x0(1) - x0(2))/(muf*Cf*Ai)==0);
u0 = double(u0);

% Box-type constraints
xmin = [0; 0];
xmax = [600; 300];
umin = 0;
umax = 0.35;
X = Polyhedron('lb', xmin, 'ub', xmax);
U = Polyhedron('lb', umin, 'ub', umax);

% Wheight matrix
Qx = 1e-4*diag([0 1]);
Rx = 1e-10;

%% Sets and set-point
Tp_sp = 109.92832; % Tf_sp = 97
rho1 = di*pi*hi*(1-exp(-Tp_sp/600))/(1-exp(-1));
syms Tf
Tf_sp = solve((-de*pi*h0/(mum*Cm*Ae)-rho1/(mum*Cm*Ae))*Tp_sp+rho1*Tf/(mum*Cm*Ae)+de*pi*nhu*Irr/(mum*Cm*Ae)+de*pi*h0*Te/(mum*Cm*Ae) == 0);
Tf_sp = double(Tf_sp);
xsp = [Tp_sp; Tf_sp];
% U
rho1 = di*pi*hi*(1-exp(-Tp_sp/600))/(1-exp(-1));
rho2 = (1-exp(-Tf_sp/300))/((1-exp(-1))*Ai);
syms u
usp = solve(-u*rho2 + rho1*(Tp_sp - Tf_sp)/(muf*Cf*Ai)==0);
usp = double(usp);

% System
A = [(-de*pi*h0/(mum*Cm*Ae)-rho1/(mum*Cm*Ae)) rho1/(mum*Cm*Ae);
        rho1/(muf*Cf*Ai) -rho1/(muf*Cf*Ai)]
B1 = [0; -rho2]

% Euler discretization method
Ad = (A*Ts) + eye(nx); B1d = B1*Ts; Cd = C; Dd = D; % Modelo en el punto inicial

% Sets
ConjX = [X.A X.b];
ConjU = [U.A U.b];

[Klqr, Plqr] = dlqr(Ad, B1d, Qx, Rx);
termSet = InvariantReg(Ad, B1d, Cd, Dd, -Klqr, ConjX, ConjU, xsp, usp);
Xterm = Polyhedron(termSet(:, 1:end-1), termSet(:, end));
Xtermx = projection(Xterm, 1:nx, 'vrep');

rho1_min=0.001;
rho1_max=di*pi*hi;
rho2_min=0.001;
rho2_max=1/Ai;

% Tp_modelo = 50; 
% rho1 = di*pi*hi*(1-exp(-Tp_modelo/600))/(1-exp(-1));
% syms Tf
% Tf_modelo = solve((-de*pi*h0/(mum*Cm*Ae)-rho1/(mum*Cm*Ae))*Tp_modelo+rho1*Tf/(mum*Cm*Ae)+de*pi*nhu*Irr/(mum*Cm*Ae)+de*pi*h0*Te/(mum*Cm*Ae) == 0);
% Tf_modelo = double(Tf_modelo);
% % U
% rho1_min = di*pi*hi*(1-exp(-Tp_modelo/600))/(1-exp(-1));
% rho2_min = (1-exp(-Tf_modelo/300))/((1-exp(-1))*Ai);

% Tp_modelo = 600; % Tf_sp = 97
% rho1 = di*pi*hi*(1-exp(-Tp_modelo/600))/(1-exp(-1));
% syms Tf
% Tf_modelo = solve((-de*pi*h0/(mum*Cm*Ae)-rho1/(mum*Cm*Ae))*Tp_modelo+rho1*Tf/(mum*Cm*Ae)+de*pi*nhu*Irr/(mum*Cm*Ae)+de*pi*h0*Te/(mum*Cm*Ae) == 0);
% Tf_modelo = double(Tf_modelo);
% % U
% rho1_max = di*pi*hi*(1-exp(-Tp_modelo/600))/(1-exp(-1));
% % rho2_max = (1-exp(-Tf_modelo/300))/((1-exp(-1))*Ai);

%% LTI vertex System 1
% Consider rho1 min y rho2 min
% Parameters
rho1 = rho1_min;
rho2 = rho2_min;

% System
A = [(-de*pi*h0/(mum*Cm*Ae)-rho1/(mum*Cm*Ae)) rho1/(mum*Cm*Ae);
        rho1/(muf*Cf*Ai) -rho1/(muf*Cf*Ai)];
B1 = [0; -rho2];

% Euler discretization method
A1d = (A*Ts) + eye(nx); B11d = B1*Ts; Cd = C; Dd = D; % Modelo en el punto inicial

[Klqr, Plqr1] = dlqr(A1d, B11d, Qx, Rx);
% termSet = InvariantReg(A1d, B11d, Cd, Dd, -Klqr, ConjX, ConjU, xsp, usp);
% Xterm = Polyhedron(termSet(:, 1:end-1), termSet(:, end));
% Xtermx = projection(Xterm, 1:nx, 'vrep');

%% LTI vertex System 2
% Consider rho1 min y rho2 max
% Parameters
rho1 = rho1_min;
rho2 = rho2_max;

% System
A = [(-de*pi*h0/(mum*Cm*Ae)-rho1/(mum*Cm*Ae)) rho1/(mum*Cm*Ae);
        rho1/(muf*Cf*Ai) -rho1/(muf*Cf*Ai)];
B1 = [0; -rho2];

% Euler discretization method
A2d = (A*Ts) + eye(nx); B12d = B1*Ts; Cd = C; Dd = D; % Modelo en el punto inicial

[Klqr, Plqr2] = dlqr(A2d, B12d, Qx, Rx);
% termSet = InvariantReg(A2d, B12d, Cd, Dd, -Klqr, ConjX, ConjU, xsp, usp);
% Xterm = Polyhedron(termSet(:, 1:end-1), termSet(:, end));
% Xtermx1 = projection(Xterm, 1:nx, 'vrep');
% Xtermx = intersect(Xtermx,Xtermx1);

%% LTI vertex System 3
% Consider rho1 max y rho2 min
% Parameters
rho1 = rho1_max;
rho2 = rho2_min;

% System
A = [(-de*pi*h0/(mum*Cm*Ae)-rho1/(mum*Cm*Ae)) rho1/(mum*Cm*Ae);
        rho1/(muf*Cf*Ai) -rho1/(muf*Cf*Ai)];
B1 = [0; -rho2];

% Euler discretization method
A3d = (A*Ts) + eye(nx); B13d = B1*Ts; Cd = C; Dd = D; % Modelo en el punto inicial

[Klqr, Plqr3] = dlqr(A3d, B13d, Qx, Rx);
% termSet = InvariantReg(A3d, B13d, Cd, Dd, -Klqr, ConjX, ConjU, xsp, usp);
% Xterm = Polyhedron(termSet(:, 1:end-1), termSet(:, end));
% Xtermx1 = projection(Xterm, 1:nx, 'vrep');
% Xtermx = intersect(Xtermx,Xtermx1);

%% LTI vertex System 4
% Consider rho1 max y rho2 max
% Parameters
rho1 = rho1_max;
rho2 = rho2_max;

% System
A = [(-de*pi*h0/(mum*Cm*Ae)-rho1/(mum*Cm*Ae)) rho1/(mum*Cm*Ae);
        rho1/(muf*Cf*Ai) -rho1/(muf*Cf*Ai)];
B1 = [0; -rho2];

% Euler discretization method
A4d = (A*Ts) + eye(nx); B14d = B1*Ts; Cd = C; Dd = D; % Modelo en el punto inicial

[Klqr, Plqr4] = dlqr(A4d, B14d, Qx, Rx);
% termSet = InvariantReg(A4d, B14d, Cd, Dd, -Klqr, ConjX, ConjU, xsp, usp);
% Xterm = Polyhedron(termSet(:, 1:end-1), termSet(:, end));
% Xtermx1 = projection(Xterm, 1:nx, 'vrep');
% Xtermx = intersect(Xtermx,Xtermx1);
asdf
figure()
plot(Xtermx, 'Color', 'g', 'Alpha', 0.05);

%% B2
B2 = [de*pi*nhu/(mum*Cm*Ae)  de*pi*h0/(mum*Cm*Ae); 0 0];
B2d = B2*Ts;

%% System Polytope
nvert=4;
A_poli=[A1d A2d A3d A4d];
B1_poli=[B11d B12d B13d B14d];
Plqr_poli = [Plqr1 Plqr2 Plqr3 Plqr4];

%% MHE
run LPVMPC_MHE

%% Run MPC controller
% Plqrmu=Plqr3;
run LPVMPC_MPC

%% Simulation
% Vector initialization for plots
Obj = []; state = x0; input = []; u_max = umax; u_min = umin; y = x0; 
W = []; Y = C*x0; X = x0; Xart = x0; mues=[]; muesMHE=[]; 
rhosreal = zeros(2,Nsim); muesreal = zeros(4,Nsim); deltasumod = zeros(1,Nsim);
Predict = zeros(2, N);         
for k = 1:N
            Predict(:, k) = [IC151(k+1, 2); TA075(k+1, 2)];
end
X_MHE = repmat(x0,1,N+1); U_MHE = repmat(u0,1,N); 
j = 1;
% mu_mhe=[0.25;0.25;0.25;0.25];
%mu_mhe=[0;0.5;0;0.5];
mu_mhe=[0;0;0;1];
disp('Iniciando...')
umpc = 0;
elapsed_time = zeros(Nsim, 1) ;  % initilaize the elapsed times 
    for i = 1:Nsim
        tk = i*Ts;
        disp(i);
        
        % Update disturbance
        Irr = IC151(i, 2);
        Te = TA075(i, 2);
        w = [Irr Te]';
        
        t_tic = tic ;   % to get time evaluated 
        
        [sol, diag] = mhe{X_MHE, U_MHE, Predict, mu_mhe}; % 
        if diag
           msg = ['Infeasible MHE at t = ', num2str(i*Ts)];
           disp(msg)
           return;
        end
        mu_mhe=sol;
        
        %mu_fuzzy = membership(x0, rho1_min, rho1_max, rho2_min, rho2_max)';
        mu_fuzzy = [1 1 1 1]';
        
        % Disturbance prediction
        %mu_mhe = [0.25;0.25;0.25;0.25];
        for k = 1:N
            Predict(:, k) = [IC151(i+k+1, 2); TA075(i+k+1, 2)];
        end
        uant = umpc;
        %[sol, diag] = mpc{x0, xsp, Predict, mu_fuzzy};
        [sol, diag] = mpc{x0, xsp, Predict, mu_mhe};
        
        if diag
            msg = ['Infeasible MPC at t = ', num2str(i*Ts)];
            disp(msg)
            return;
        end
        umpc = sol{1}; 
        obj = sol{2};
        
        t_tic = toc(t_tic) ;              % get time elapsed
        elapsed_time(i) = t_tic ;   % store the time elapsed for the run

        mues = [mues mu_mhe];
        mues = [mues mu_fuzzy];
        muesMHE =[muesMHE mu_mhe];
        rhosreal(1,i) = di*pi*hi*(1-exp(-X(1,end)/600))/(1-exp(-1));
        rhosreal(2,i) = (1-exp(-X(2,end)/300))/(Ai*(1-exp(-1)));
        
        muesreal(1,i) = ((rho1_max-rhosreal(1,i))/(rho1_max-rho1_min))*((rho2_max-rhosreal(2,i))/(rho2_max-rho2_min));
        muesreal(2,i) = ((rho1_max-rhosreal(1,i))/(rho1_max-rho1_min))*((rhosreal(2,i)-rho2_min)/(rho2_max-rho2_min));
        muesreal(3,i) = ((rhosreal(1,i)-rho1_min)/(rho1_max-rho1_min))*((rho2_max-rhosreal(2,i))/(rho2_max-rho2_min));
        muesreal(4,i) = ((rhosreal(1,i)-rho1_min)/(rho1_max-rho1_min))*((rhosreal(2,i)-rho2_min)/(rho2_max-rho2_min));
        

        
        % Continuous-time simulation (reality)
        [tsim, x] = ode45(@ST, [0 Ts], x0, options, umpc, w);
        x0 = x(end, :)';                 % Next state

        X = [X x(2:end, :)'];          % State vector
        Y = [Y C*x(2:end, :)'];       % Output vector
        t = [t tk+tsim(2:end)'];    % Time vector
        input = [input; umpc];     % Input vector
        Obj = [Obj obj];                % Objective vector
        W = [W w];                        % Disturbance vector
        
        % Update MHE
        X_MHE = [X_MHE(:,2:end) x0];
        U_MHE = [U_MHE(:,2:end) umpc];
        
        deltasumod(1,i) = abs(umpc-uant); 
        
    end
    
    
figure;
subplot(3,1,1);
plot(Ts:Ts:Ts*Nsim,0.61754*muesreal(1,:)/0.411435); hold on; grid on;
plot(Ts:Ts:Ts*Nsim,muesMHE(1,:));
plot(Ts:Ts:Ts*Nsim,0.11447*muesreal(2,:)/0.32);
plot(Ts:Ts:Ts*Nsim,muesMHE(2,:));
plot(Ts:Ts:Ts*Nsim,0.12216*muesreal(3,:)/0.1517);
plot(Ts:Ts:Ts*Nsim,muesMHE(3,:));
plot(Ts:Ts:Ts*Nsim,0.146*muesreal(4,:)/0.118);
plot(Ts:Ts:Ts*Nsim,muesMHE(4,:));
subplot(3,1,2);
plot(Ts:Ts:Ts*Nsim,ones(1,Nsim)*rho1_min); hold on; grid on;
plot(Ts:Ts:Ts*Nsim,rhosreal(1,:));
plot(Ts:Ts:Ts*Nsim,ones(1,Nsim)*rho1_max);
subplot(3,1,3);
plot(Ts:Ts:Ts*Nsim,ones(1,Nsim)*rho2_min); hold on; grid on;
plot(Ts:Ts:Ts*Nsim,rhosreal(2,:));
plot(Ts:Ts:Ts*Nsim,ones(1,Nsim)*rho2_max);

% Perfonmance indices
error = abs(X(2, :) - xsp(2));
IAE = trapz(t, abs(error));
ISE = trapz(t, error.^2);
ITAE = trapz(t, t.*abs(error));
msg = ['IAE = ', num2str(IAE)];
disp(msg)
msg = ['TV = ', num2str(sum(deltasumod(1,50:end)))];
disp(msg)
msg = ['ITAE = ', num2str(ITAE)];
disp(msg)
msg = ['ISE = ', num2str(ISE)];
disp(msg)

error = abs(X(2, 1:7000) - xsp(2));
IAEtrack = sum(error)/length(error);
error = abs(X(2, 7001:end) - xsp(2));
IAEpert = sum(error)/length(error);

time_avg = mean(elapsed_time) ;
msg = ['Mean time = ', num2str(time_avg)];
disp(msg)
time_avg = max(elapsed_time) ;
msg = ['Max time = ', num2str(time_avg)];
disp(msg)
time_avg = min(elapsed_time) ;
msg = ['Min time = ', num2str(time_avg)];
disp(msg)

save run.mat

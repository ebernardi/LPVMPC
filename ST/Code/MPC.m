%% Forward MPC QP
yalmip('clear')

%% Controller definition
x = sdpvar(nx*ones(1, N_MPC+1), ones(1, N_MPC+1));
u = sdpvar(nu*ones(1, N_MPC), ones(1, N_MPC));
us = sdpvar(nu, 1);
xs = sdpvar(nx, 1);
mu = sdpvar(M, 1);
Pred = sdpvar(nx, N_MPC);

% Artificial variables
xa = sdpvar(nx, 1); 
ua = sdpvar(nu, 1);

objective = 0; constraints = [];
gamma = 1e2*Plqr;

A = [sys.Ad]*kron(mu, eye(nx));
B = [sys.Bd]*kron(mu, eye(nu));
Bw = [sys.Bwd]*kron(mu, eye(nd));
delta = [sys.deltad]*mu;

% Stage constraints and objective
for k = 1:N_MPC
    % Objective
    % With auxiliary variable
    objective = objective + (x{k}-xa)'*Qx*(x{k}-xa); 
    objective = objective + (u{k}-ua)'*Ru*(u{k}-ua);
% 	% Without auxiliary variable
%     objective = objective + (x{k}-xs)'*Qx*(x{k}-xs); 
%     objective = objective + (u{k}-us)'*Ru*(u{k}-us);

    % Dynamic constraint    
    constraints = [constraints, x{k+1} == A*x{k} + B*u{k} + Bw*Pred(:, k) + delta];

    % Box-type constraint
    constraints = [constraints, umin <= u{k} <= umax];
    constraints = [constraints, xmin <= x{k} <= xmax];
end

%% BMI 
% With auxiliary variable
% Terminal cost
objective = objective + (xa-xs)'*gamma*(xa-xs);
objective = objective + (x{N_MPC+1}-xa)'*Pbmi*(x{N_MPC+1}-xa); 

% Terminal constraint
constraints = [constraints, xa == A*xa + B*ua + Bw*Pred(:, N_MPC) + delta]; % Artificial variables equilibirum condition
constraints = [constraints, (x{N_MPC+1}-xa)'*Wbmi*(x{N_MPC+1}-xa) <= 1];

% % Without auxiliary variable
% % Terminal cost
% objective = objective + (x{N_MPC+1}-xs)'*Pbmi*(x{N_MPC+1}-xs);
% % Terminal constraint
% constraints = [constraints, (x{N_MPC+1}-xs)'*Wbmi*(x{N_MPC+1}-xs) <= 1];

% Defining the parameters in, and the solution
parameters = {x{1}, xs, us, Pred, mu};
solution = {u{1}, [x{:}], objective};

% Options for Optimizer  
options = sdpsettings('solver', 'quadprog');
options.verbose = 0;                             % 0 to none, 2 to debug, 2+ for more

mpc = optimizer(constraints, objective, options, parameters, solution);
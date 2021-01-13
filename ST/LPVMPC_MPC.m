%% Forward MPC QP
yalmip('clear')

%% Controller definition
x = sdpvar(nx*ones(1, N+1), ones(1, N+1));
u = sdpvar(nu*ones(1, N), ones(1, N));
xs = sdpvar(nx, 1);
mu = sdpvar(nvert,1);
Pred = sdpvar(nx, N);

% Artificial variables
xa = sdpvar(nx, 1); 
ua = sdpvar(nu, 1);

objective = 0; constraints = [];
% gamma = 10*Plqr;

Amu = A_poli*kron(mu,eye(nx));            % A(mu) = A1*mu1 + A2*mu2 + ...
Bmu = B1_poli*kron(mu,eye(nu));            % B(mu) = B1*mu1 + B2*mu2 + ...
Plqrmu = Plqr_poli*kron(mu,eye(nx));            % Plqr(mu) = Plqr*mu1 + Plqr*mu2 + ...

% Stage constraints and objective
for k = 1:N
    % Objective
    objective = objective + (x{k}-xa)'*Qx*(x{k}-xa); 
    objective = objective + (u{k}-ua)'*Rx*(u{k}-ua);

    % Dynamic constraint    
    constraints = [constraints, x{k+1} == Amu*x{k} + Bmu*(u{k}) + B2d*Pred(:, k)];

    % Box-type constraint
    constraints = [constraints, umin <= u{k} <= umax];
    constraints = [constraints, xmin <= x{k} <= xmax];
end

gamma = 10*Plqrmu;
% Terminal constraints
objective = objective + (x{N+1}-xa)'*Plqrmu*(x{N+1}-xa); 
objective = objective + (xa-xs)'*gamma*(xa-xs);                % Terminal cost

constraints = [constraints, xa == Amu*xa + Bmu*ua + B2d*Pred(:, N)]; % Artificial variables equilibirum condition
% constraints = [constraints, Xtermx.A*x{N+1}<=Xtermx.b];

% Defining the parameters in, and the solution
parameters = {x{1}, xs, Pred, mu};
solution = {u{1}, objective};

% Options for Optimizer  
options = sdpsettings('solver', 'quadprog');
options.verbose = 0;                             % 0 to none, 2 to debug, 2+ for more

mpc = optimizer(constraints, objective, options, parameters, solution);
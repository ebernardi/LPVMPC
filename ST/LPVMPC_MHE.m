%% Backward MHE QP
yalmip('clear')

%% MHE definition
x    = sdpvar(nx*ones(1, N+1), ones(1, N+1));
xru  = sdpvar(nx*ones(1, N+1), ones(1, N+1));
nux  = sdpvar(nx*ones(1, N+1), ones(1, N+1));
u    = sdpvar(nu*ones(1, N), ones(1, N));
mu   = sdpvar(nvert,1);
muan = sdpvar(nvert,1);
numu = sdpvar(nvert,1);
e    = sdpvar(nx*ones(1, N+1), ones(1, N+1));
Pred = sdpvar(nx, N);

objective = 0; constraints = [];

Amu = A_poli*kron(mu,eye(nx));            % A(mu) = A1*mu1 + A2*mu2 + ...
Bmu = B1_poli*kron(mu,eye(nu));            % B(mu) = B1*mu1 + B2*mu2 + ...

Qe    = 3*52500;
Qnux  = 10;
Qnumu = 5*10000;


% Stage constraints and objective
for k = 1:N
    % Objective
    objective = objective + (e{k})'*Qe*(e{k}) + (nux{k})'*Qnux*(nux{k}); 

    % Dynamic constraint    
    constraints = [constraints, e{k+1} == xru{k+1} - (Amu*x{k} + Bmu*u{k} + B2d*Pred(:, k))];
    constraints = [constraints, xru{k} == x{k} + nux{k}];   
end

constraints = [constraints, mu == muan + numu];
constraints = [constraints, sum(mu) == 1];
constraints = [constraints, zeros(nvert,1) <= mu <=ones(nvert,1)];

objective = objective + (numu)'*Qnumu*(numu);

% Defining the parameters in, and the solution
parameters = {[x{:}], [u{:}], Pred, [muan]}; % 
solution = {mu};

% Options for Optimizer  
options = sdpsettings('solver', 'gurobi');
options.verbose = 0;                             % 0 to none, 2 to debug, 2+ for more

mhe = optimizer(constraints, objective, options, parameters, solution);
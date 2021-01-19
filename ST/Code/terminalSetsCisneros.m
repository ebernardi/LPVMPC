%% Terminal ingredients [Cisneros et al.] (2016)
yalmip('clear')

% Define variables
Y = sdpvar(nx, nx, 'symmetric');
X = sdpvar(nu, nx, 'full');
t_opt = sdpvar(1);
alpha_tilde = sdpvar(1);

% Define constraints 
constraints = [];

constraints = [constraints, Y >= 0];
for i = 1:M
    constraints = [constraints, [Y (Y*sys(i).Ad'+X'*sys(i).Bd') Y X'; ...
                                                (sys(i).Ad*Y+sys(i).Bd*X) Y zeros(nx, nx) zeros(nx, nu); ...
                                                Y zeros(nx, nx) Qx^(-1) zeros(nx, nu); ...
                                                X zeros(nu, nx) zeros(nu, nx) Ru^(-1)] >= 0];
end
% Define an objective
objective = [];

% Solve the problem
sol = optimize(constraints, objective);

Pbmi = inv(value(Y));
K = value(X)/value(Y);

% Second problem
yalmip('clear')

% Define variables
Z = sdpvar(nx, nx, 'symmetric');

% Define constraints 
constraints = [];
constraints = [constraints, Z >= 0];
for i = 1:M
    constraints = [constraints, [-Z Z*(sys(i).Ad+sys(i).Bd*K)'; ...
                                                (sys(i).Ad+sys(i).Bd*K)*Z -Z] <= 0];
end
constraints = [constraints, K*Z*K' <= u_abs^2];

% Define an objective
objective = -logdet(Pbmi);

% Solve the problem
sol = optimize(constraints, objective);

Wbmi = inv(value(Z));

% ellipse(Wbmi, Xsp, 20, 'black', '-')
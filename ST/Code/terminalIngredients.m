%% [Kothare] (1996)
% First LMI problem
yalmip('clear')

Y = sdpvar(nx, nx, 'symmetric');
L = sdpvar(nu, nx, 'full');
sigma = sdpvar(1);

% Define constraints 
constraints = [];

constraints = [constraints, Y >= 0];
constraints = [constraints, sigma >= 0];
for i = 1:M
    constraints = [constraints, [Y Y*sys(i).Ad'+L'*sys(i).Bd' Y*Qx^(1/2) L'*Ru^(1/2); ...
                                                sys(i).Ad*Y+sys(i).Bd*L Y zeros(nx, nx) zeros(nx, nu); ...
                                                Qx^(1/2)*Y zeros(nx, nx) sigma*eye(nx, nx) zeros(nx, nu); ...
                                                Ru^(1/2)*L zeros(nu, nx) zeros(nu, nx) sigma*eye(nu, nu)] >= 0];
end
% Define an objective
objective = sigma;

% Solve the problem
ops = sdpsettings('solver', 'sedumi');
sol = optimize(constraints, objective, ops);

string = yalmiperror(sol.problem);
msg = ['Output => ', string];
disp(msg)

if sol.problem
	return;
end

sigma = value(sigma);
Pbmi = eye(size(value(Y)))*value(sigma)/value(Y);
K = value(L)/value(Y);

%% [Gonzales] (2016)
% Second LMI problem
yalmip('clear')

% Define variables
Z = sdpvar(nx, nx, 'symmetric');

% Define constraints 
constraints = [];
constraints = [constraints, Z >= 0];
for i = 1:M
    constraints = [constraints, [Z (sys(i).Ad+sys(i).Bd*K)*Z; ...
                                                Z*(sys(i).Ad+sys(i).Bd*K)' Z] >= 0];
end
for j = 1:nu
        constraints = [constraints, K(j, :)*Z*K(j, :)' >= (u_abs(j, :)^2)];
end

% Define an objective
objective = -geomean(Z);

ops = sdpsettings('solver', 'sedumi');
sol = optimize(constraints, objective, ops);

Wbmi = eye(size(value(Z)))/value(Z);

figure(1)
hold on
ellipse(Wbmi, Xsp, 20, 'red', '-')
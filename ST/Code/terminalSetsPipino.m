%% Terminal Ingredients
yalmip('clear')

%% Terminal cost weight 
Y = sdpvar(nx, nx, 'symmetric');
L = sdpvar(nu, nx, 'full');
sigma = sdpvar(1);
constraints = [];

constraints = [constraints, Y >= 0];
for i = 1:M
    constraints = [constraints, [Y Y*sys(i).Ad'+L'*sys(i).Bd' Y*Qx^(1/2) L'*Ru^(1/2); ...
                                                sys(i).Ad*Y+sys(i).Bd*L Y zeros(nx, nx) zeros(nx, nu); ...
                                                Qx^(1/2)*Y zeros(nx, nx) sigma*eye(nx, nx) zeros(nx, nu); ...
                                                Ru^(1/2)*L zeros(nu, nx) zeros(nu, nx) sigma*eye(nu, nu)]];
end
optimize(constraints, sigma);
Plqr = eye(size(value(Y)))/(value(Y)/value(sigma));
Klqr = value(L)/value(Y)

%% Terminal set constraint
Z = Polyhedron('lb', [xmin; umin], 'ub', [xmax; umax]); % Extended set

X = projection(Z, 1:nx); X = minHRep(X);
U = projection(Z, nx+1:nx+nu); U = minHRep(U);
setX = [X.A X.b];
setU = [U.A U.b];

for i = 1:M
    sys(i).Plqr = Plqr;
    sys(i).Klqr = (Ru+sys(i).Bd'*Plqr*sys(i).Bd)\(sys(i).Bd'*Plqr*sys(i).Ad);
    
    Oinf = InvariantReg(sys(i).Ad, sys(i).Bd, sys(i).Cd, sys(i).Dd, -sys(i).Klqr, setX, setU, xsp, usp);
    invariant_set = Polyhedron(Oinf(:, 1:end-1), Oinf(:, end));
    invariant_set = minHRep(invariant_set);
    if i == 1
        Xf = invariant_set;
    else
        Xf = intersect(Xf, invariant_set);
    end
end
function Oinft = InvariantReg(A, B, C, D, Klqr, Xo, Uo, xs, us)

[nx, nu] = size(B);
ny = length(C(:, 1));

if isempty(Xo)
   nrx = 0;
else
   nrx = length(Xo(:, 1));
end

if isempty(Uo)
   nru = 0;
else
   nru = length(Uo(:, 1));
end

% Generamos un conjunto con centro en el origen
X = [Xo(:, 1:end-1), Xo(:, end)-Xo(:, 1:end-1)*xs]; % Conjunto de restricciones de estados
U = [Uo(:, 1:end-1), Uo(:, end)-Uo(:, 1:end-1)*us]; % Conjunto de restricciones de manipuladas

% display('Conjunto Terminal');
Acl = A+B*Klqr; %Matriz A del sistema extendido de lazo cerrado

% Conjunto de restricciones del sistema extendido
W = zeros(nrx+nru, nx+1); 
W(1:nrx, :) = [X(:, 1:nx), X(:, nx+1)];
W(nrx+1:nrx+nru, :) = [U(:, 1:nu)*Klqr, U(:, nu+1)];

Oi = W;
maxiteracion = 101;

for i=2:maxiteracion
    Oir = Oi;
    Oi = remred([Oir; [W(:, 1:end-1)*Acl^(i-1), W(:, end)]]);
    if issubset(Oir, Oi)
%          fprintf('Invariante con numero de pasos = %d \n', i);
        break;
    end
end

Oinf = Oi;

% Lo llevamos nuevamente a xs y us
Oinft = [Oinf(:, 1:end-1), Oinf(:, end)+Oinf(:, 1:end-1)*xs];
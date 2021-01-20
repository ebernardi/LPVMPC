function fval = issubset(U, V, tol)
% ISSUBSET True if polyhedron U is a subset of or equal to polyhedron V
%
% Usage: 
% fval = issubset(U, V, tol, display)

P = polytope(U(:, 1:end-1), U(:,end));
Q = polytope(V(:, 1:end-1), V(:,end));

if nargin == 2
    flag = le(P, Q);
elseif nargin == 3
    % Amplio V en tol para ver si cabe
    Options.abs_tol = tol;
    Options.rel_tol = tol;    
    flag = le(P, Q, Options);
end

if flag
    fval = 1;
else
    fval = 0;
end
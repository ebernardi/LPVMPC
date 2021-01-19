function [Anew, bnew] = remred(A, b)
% Funcion de kerrigan a mpt

if nargin == 1
    P = A;
    A = P(:, 1:end-1);
    b = P(:, end);
end

P = Polyhedron(A, b);
PR = minHRep(P);
   
Anew = PR.A;
bnew = PR.b;

if nargout == 1
    Anew = [Anew, bnew];
end
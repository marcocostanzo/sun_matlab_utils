function [isPD,isSPD] = isPositiveDefinite(A, tolerance)
% check positive definite
% isPD = is Positive Definite
% isSPD = is semiPositive definite

    d = eig(A);
    
    if nargin < 2 || isempty(tolerance)
        tolerance = length(d)*eps(max(d));
    end
    
    isSPD = all(d) > -tolerance;
    isPD = all(d) > tolerance;

end

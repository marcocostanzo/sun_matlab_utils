function issym = isSymmetric(A,tolerance)
%ISSYMMETRIC Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2 || isempty(tolerance)
	tolerance = length(A)*eps(max(A,[],'all'));
end

issym = norm(A-A') < tolerance;

end


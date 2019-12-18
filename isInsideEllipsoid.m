function inside = isInsideEllipsoid(x,A,C)
%ISINSIDEELLIPSOID Check if x is inside the ellipsoid defined by A, and C
%   (x-C)'*A*(x-C) <= 1
%   Detailed explanation goes here

inside = (x-C)'*A*(x-C) - 1 <= 0;

end


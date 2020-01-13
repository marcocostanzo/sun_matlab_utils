function [d, x, y, iters] = ellipsoid_distance(A1,center1,A2,center2,c1,c2,epsilon,max_iter, do_checks)
%ELLIPSOID_DISTANCE Compute the distance between two ellipsoids
%   the ellipsoids are in the form x'*(A-center)*x < 1
%
%
% Algorithm from:
% Lin, A., & Han, S.-P. (2002). 
% On the Distance between Two Ellipsoids. 
% SIAM Journal on Optimization, 13(1), 298–308. 
% doi:10.1137/s1052623401396510 
%

%% Assertions

if nargin < 9
    do_checks = true;
end

if do_checks

assert( ...
    all(size(A1) == size(A2)) ...
    && ...
    size(A1,1) == size(A1,2) ...
    && ...
    size(A1,1) == length(center1) ...
    , 'ellipsoid_distance:size_mismatch', 'size_mismatch')

%% Check h.
% check if matrix A1 A2 are symmetric positive definite
if ~isSymmetric(A1)
    error('ellipsoid_distance:NOT_Symmetric','Matrix A1 is not Symmetric')
end
[isPD,isSPD] = isPositiveDefinite(A1);
if ~isPD
    if isSPD % check semipositive definite to give a warning
        warning('ellipsoid_distance:SPD','Matrix A1 is SPD')
    else
        error('ellipsoid_distance:NOT_SPD','Matrix A1 is not SPD')
    end
end
if ~isSymmetric(A2)
    error('ellipsoid_distance:NOT_Symmetric','Matrix A2 is not Symmetric')
end
[isPD,isSPD] = isPositiveDefinite(A2);
if ~isPD
    if ~isSPD % check semipositive definite to give a warning
        warning('ellipsoid_distance:SPD','Matrix A2 is SPD')
    else
        error('ellipsoid_distance:NOT_SPD','Matrix A2 is not SPD')
    end
end

%% Compute the paper form 0.5*x'*A*x + b'*x + a <= 0
% Not needed in this version
% a1_bar = center1'*A1*center1;
% b1_bar = -2*A1*center1;
% A1_bar = 2*A1;
% 
% a2_bar = center2'*A2*center2;
% b2_bar = -2*A2*center2;
% A2_bar = 2*A2;

%% Init
if nargin < 8 || isempty(max_iter)
    max_iter = 3000;
end
if nargin < 7 || isempty(epsilon)
    epsilon = deg2rad(0.1);
end
%Initial points are c1€E1 and C2€E2
%If not provided the centers be used
if nargin < 6 || isempty(c2)
    c2 = center2; %-A2\b2;
end
if nargin < 5 || isempty(c1)
    c1 = center1; %-A1\b1;
end

%% Check initial value of c1 and c2
if ~isInsideEllipsoid(c1,A1,center1)
    error('ellipsoid_distance:c_not_in_E','Initial value of c1 is outside E')
end
if ~isInsideEllipsoid(c2,A2,center2)
    error('ellipsoid_distance:c_not_in_E','Initial value of c2 is outside E')
end

end %end do checks

%% iters
iters = 1;

b_stop = false;

while(~b_stop)
    
    % solve eq (2.2) of paper to get t1 and t2
    t1 = getStepSize(c1,c2,A1,center1,true);
    t2 = getStepSize(c1,c2,A2,center2,false);
    
    if t2 <= t1
        % DONE! distance is zero
        d = 0;
        x = nan(size(A1,1),1);
        y = x;
        return;
    end
    
    % update x and y
    x = c1 + t1*(c2-c1);
    y = c1 + t2*(c2-c1);
    
    % compute angle theta1 and theta2
    % NB in paper is:
    % theta1 = theta_fun(y-x,A1_bar*x+b1_bar);
    % theta2 = theta_fun(x-y,A2_bar*y+b2_bar);
    % BUT:
    % A_bar*x+b = 2*A*(x-C)
    % I define vec_i = A_bar_i*[x||y]+b_i = 2*A_i*(x-C_i)
    vec_1 = 2*A1*(x-center1);
    vec_2 = 2*A2*(y-center2);
    theta1 = theta_fun(y-x,vec_1);
    theta2 = theta_fun(x-y,vec_2);
    
    if (theta1<epsilon) && (theta2<epsilon)
        % DONE
        d = norm(x-y);
        return;
    end
    
    % Update c1 c2 for the next iteration
     %gamma1 and gamma2 are chosen as (2.1)
    gamma1 = 1/norm(2*A1);
    gamma2 = 1/norm(2*A2);
    c1 = x - gamma1*(vec_1);
    c2 = y - gamma2*(vec_2);
    
    % update number of iter
    iters = iters + 1;
    
    % update b_stop
    b_stop = iters > max_iter;
    
end

% max iter reached
d = norm(x-y);

end


function theta = theta_fun(x,y)

    theta = acos( x'*y / ( norm(x)*norm(y) ) );

end

function t = getStepSize(c1,c2,A,center,use_max)

    coeff = zeros(1,3);
    
    coeff(1) = (c2-c1)'*A*(c2-c1);
    coeff(2) = 2*(c2-c1)'*A*(c1-center);
    coeff(3) = (c1-center)'*A*(c1-center) - 1;
    
%     %As in paper
%     coeff(1) = 0.5*(c2-c1)'*A*(c2-c1);
%     coeff(2) = (A*c1 + b)'*(c2-c1);
%     coeff(3) = (c1'*A + b')*c1 + a;

    D = coeff(2)^2 - 4*coeff(1)*coeff(3);
    
    r = [
        (-coeff(2) - sqrt(D))/(2*coeff(1));
        (-coeff(2) + sqrt(D))/(2*coeff(1))
    ];
%     r = roots(coeff);
    
    r = r(imag(r)==0 & r>=0 & r<=1);

    if use_max 
        if isempty(r)
            t = 1;
        else
            t = r(1);
        end
    else
        if isempty(r)
            t = 0;
        else
            t = r(1);
        end
    end

end






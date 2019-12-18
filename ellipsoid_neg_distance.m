function [d, x, y] = ellipsoid_neg_distance(A1,center1,A2,center2,c1,c2,epsilon,max_iter)
%ELLIPSOID_DISTANCE Compute the distance between two ellipsoids
%   Detailed explanation goes here
%
%
% Algorithm from:
% Lin, A., & Han, S.-P. (2002). 
% On the Distance between Two Ellipsoids. 
% SIAM Journal on Optimization, 13(1), 298�308. 
% doi:10.1137/s1052623401396510 
%

%% Assertions

assert( ...
    all(size(A1) == size(A2)) ...
    && ...
    size(A1,1) == size(A1,2) ...
    && ...
    size(A1,1) == length(center1) ...
    , 'ellipsoid_distance:size_mismatch')

%% Check h.
% check if matrix A1 A2 are symmetric positive definite
if ~fastcheck_SymPD(A1)
    if ~check_SPD(A1) % check semipositive definite to give a warning
        warning('ellipsoid_distance:SPD','Matrix A1 is SPD')
    else
        error('ellipsoid_distance:NOT_SPD','Matrix A1 is not SPD')
    end
end
if ~fastcheck_SymPD(A2)
    if ~check_SPD(A2) % check semipositive definite to give a warning
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
    max_iter = 1000;
end
if nargin < 7 || isempty(epsilon)
    epsilon = deg2rad(0.1);
end
%Initial points are c1�E1 and C2�E2
%If not provided the centers be used
if nargin < 6 || isempty(c2)
    c2 = center2; %-A2\b2;
end
if nargin < 6 || isempty(c1)
    c1 = center1; %-A1\b1;
end

%% Check initial value of c1 and c2
if ~isInsideEllipsoid(c1,A1,center1)
    error('ellipsoid_distance:c_not_in_E','Initial value of c1 is outside E')
end
if ~isInsideEllipsoid(c2,A2,center2)
    error('ellipsoid_distance:c_not_in_E','Initial value of c2 is outside E')
end

%% Check hip
% c1�E1 and c2�E2

% iters
iter = 1;

b_stop = false;

while(~b_stop)
    
    % update x,y as max distance on the segment
    [x,y] = getMaxDistancePoints(c1,c2,A1,center1,A2,center2);
    
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
        break;
    end
    
    % Update c1 c2 for the next iteration
     %gamma1 and gamma2 are chosen as (2.1)
    gamma1 = 1/norm(2*A1);
    gamma2 = 1/norm(2*A2);
    c1 = x - gamma1*(vec_1);
    c2 = y - gamma2*(vec_2);
    
    % update number of iter
    iter = iter + 1;
    
    % update b_stop
    b_stop = iter > max_iter;
    
end

d = norm(x-y);
iter

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
    
    r = roots(coeff);
    r = r(r>=0 & r<=1);

    if use_max 
        if isempty(r)
            t = 1;
        else
            t = r;
        end
    else
        if isempty(r)
            t = 0;
        else
            t = r;
        end
    end

end

function [x,y] = getMaxDistancePoints(c1,c2,A1,center1,A2,center2)
    coeff = zeros(1,3);
    
    coeff(1) = (c2-c1)'*A1*(c2-c1);
    coeff(2) = 2*(c2-c1)'*A1*(c1-center1);
    coeff(3) = (c1-center1)'*A1*(c1-center1) - 1;
    r1 = roots(coeff);
    
    coeff(1) = (c2-c1)'*A2*(c2-c1);
    coeff(2) = 2*(c2-c1)'*A2*(c1-center2);
    coeff(3) = (c1-center2)'*A2*(c1-center2) - 1;
    r2 = roots(coeff);
    
    Points = [ c1 + r1(1)*(c2-c1), c1 + r1(2)*(c2-c1), ...
               c1 + r2(1)*(c2-c1), c1 + r2(2)*(c2-c1)];
           
    Feasible_Points = [];
    for i=1:size(Points,2)  
        %check if feasible
        if (isInsideEllipsoid(Points(:,i),A1,center1) ...
         && isInsideEllipsoid(Points(:,i),A2,center2))
            Feasible_Points = [Feasible_Points Points(:,i)];
        end
    end
    
    best
    
    X = zeros(size(A1,1),4);
    Y = X;
    
    X(:,1) = c1 + r1(1)*(c2-c1);
    X(:,2) = c1 + r1(2)*(c2-c1);
    
    Y(:,1) = c1 + r2(1)*(c2-c1);
    Y(:,3) = c1 + r2(2)*(c2-c1);
    
    % use only feasible solutions
    if ~(isInsideEllipsoid(X(:,1),A1,center1) ...
         && isInsideEllipsoid(X(:,1),A2,center2))
        X(:,1) = nan;
    end
    if ~(isInsideEllipsoid(X(:,2),A1,center1) ...
         && isInsideEllipsoid(X(:,2),A2,center2))
        X(:,2) = nan;
    end
    if ~(isInsideEllipsoid(Y(:,1),A1,center1) ...
         && isInsideEllipsoid(Y(:,1),A2,center2))
        Y(:,1) = nan;
    end
    if ~(isInsideEllipsoid(Y(:,3),A1,center1) ...
         && isInsideEllipsoid(Y(:,3),A2,center2))
        Y(:,3) = nan;
    end
    
    % Duplicate for comparison
    X(:,3) = X(:,1);
    X(:,4) = X(:,2);
    
    Y(:,2) = Y(:,1);
    Y(:,4) = Y(:,3);
    
    D = vecnorm(X-Y);
    
    [~,index] = max(D);
    
    x = X(:,index);
    y = Y(:,index);    
    
end

function isSymPD = fastcheck_SymPD(A)

    try chol(A)
        isSymPD = true;
    catch ME
        isSymPD = false;
    end

end

function isSPD = check_SPD(A, tol)

    d = eig(A);
    
    if nargin < 2 || isempty(tol)
        tol = length(d)*eps(max(d));
    end
    
    isSPD = all(d) > -tol;

end


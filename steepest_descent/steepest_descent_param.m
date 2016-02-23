function [x,p] = steepest_descent_param(f, x00, p0, x_obj, tol, kmax, df, hxx, hpx)
%STEEPEST_DESCENT_PARAM inverse method of steepest descent
%   [x,p]
%   (f, x00, p0, x_obj, tol, kmax, df, hxx, hpx)
%   f : objective function to minimize
%   x0 : initial paramters (arguments for the function)
%   tol : paramter for convergence (condition < tolerance)
%   kmax : maximum iteration time
%   df : cell of derivatives (df/dx1, df/dx2, df/dx3, ...)
%   hxx : hessian, cell of cell
%   hpx : hessian, cell of cell

%%
num = true;
if (nargin == 9)
    num = false;
end

%%
n = length(x_obj);
m = length(p0);
fx0 = @(x) f(x, p0);
x0 = steepest_descent(fx0, x00, tol, kmax); % starting point for the new min problem
xk = x0;
pk = p0;
err = 100;
k = 0;
alpha0 = 100;
alphak = alpha0;
grad_old = zeros(size(pk));
grad = grad_old;

%%
while (k < kmax && err > tol)
    xk_old = xk;
    pk_old = pk;
    if num
        xpk = [xk_old, pk_old];
        fxp = @(x) f(x(1:n), x(n+1:n+m));
        hess = hessian(fxp,xpk);
        Dxx = hess(1:n,1:n);
        Dpx = hess(1:n,n+1:n+m);
        grad_x_p = - Dxx \ Dpx;
    else
        Dxx = hxx(xk_old, pk_old);
        Dpx = hpx(xk_old, pk_old);
        %{
        Dxx = zeros(n, n);
        Dpx = zeros(n, m);
        for i = 1:n
            for j = 1:n
                Dxx(i, j) = hxx{i, j}(xk_old, pk_old);
            end
            for j = 1:m
                Dpx(i, j) = hpx{i, j}(xk_old, pk_old);
            end
        end
        %}
        grad_x_p = - Dxx \ Dpx;
    end
    for i = 1:m
        x = 0;
        for j = 1:n
            x = x + 2 * (xk(j) - x_obj(j)) * grad_x_p(j, i);
        end
        grad(i) = x;
    end
    dir = -grad;

    if num
        f_p = @(p) objective_function(f, p, xk_old, x_obj, tol / 10000, kmax);
    else
        f_p = @(p) objective_function(f, p, xk_old, x_obj, tol / 10000, kmax, df);
    end
    alphak = backtracking_ls(f_p, pk_old, dir', grad', alphak*2, 0.5, 0.1, kmax); % create step 
    disp(alphak);
    %alphak = 0.1;
    pk = pk + alphak * dir;
    disp(pk)
    fx = @(x) f(x,pk); % mis  jour energie
    xk = steepest_descent(fx, xk, tol, kmax);
        
    % CONVERGENCE
    % err = cc_grad(grad); % condition of gradient
    % err = cc_func(f(xk_old), f(xk)); % condition of objective function difference
    % err = cc_dgrad(grad_old, grad); % condition of gradient difference
    err = cc_var(xk_old, xk); % condition of parameters difference
    % err = min(cc_func(f(xk_old), f(xk)), min(cc_grad(grad), cc_var(xk_old, xk))); % choose value for convergence
    
    % COUNTUP
    k = k + 1; % increase count
            
end

%%
x = xk;
p = pk;

end
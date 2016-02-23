function [x] = steepest_descent(f, x0, tol, kmax, df)
%STEEPEST_DESCENT method of steepest descent
%   [x]
%   (f, x0, tol, kmax, df)
%   f : objective function to minimize
%   x0 : initial paramters (arguments for the function)
%   tol : paramter for convergence (condition < tolerance)
%   kmax : maximum iteration time
%   df : cell of derivatives (df/dx1, df/dx2, df/dx3, ...)

%%
num = true;
if (nargin == 5)
    num = false;
end

%%
xk = x0; % initializing
err = 10000; % err for convergence, initializing
k = 0; % counter
alphak = 100;  % initializing
[grad_old, ~] = gradest(f,xk);  % initializing
if not(num)
    for i = 1:length(x0)
        grad_old(i) = df{i}(xk);  % initializing
    end
end
grad = grad_old;  % initializing

%%
while (k < kmax && err > tol)
    
    % UPDATE
    xk_old = xk; % record former position
    grad_old = grad; % record former gradient 
    pk = - grad_old'; % search direction 
    alphak = backtracking_ls(f, xk, pk, grad_old', alphak*2, 0.5, 0.5, 1000); % create step 
    xk = xk_old + alphak*pk'; % update position
    if num 
        [grad, ~] = gradest(f,xk); % update gradient
    else
        for i = 1:length(x0)
            grad(i) = df{i}(xk); % update gradient
        end
    end
    
    % CONVERGENCE
    err = cc_grad(grad); % condition of gradient
    % err = cc_func(f(xk_old), f(xk)); % condition of objective function difference
    % err = cc_dgrad(grad_old, grad); % condition of gradient difference
    % err = cc_var(xk_old, xk)); % condition of parameters difference
    % err = min(cc_func(f(xk_old), f(xk)), min(cc_grad(grad), cc_var(xk_old, xk))); % choose value for convergence
    
    % COUNTUP
    k = k + 1; % increase count
    
end

%%
x = xk; % return value

end


function [x, x_list, cc_grad_list, cc_func_list, cc_dgrad_list, cc_var_list, alpha_list] = steepest_descent_record(f, x0, tol, kmax, df)
%STEEPEST_DESCENT_RECORD method of steepest descent
%   [x, x_list, cc_grad_list, cc_func_list, cc_dgrad_list, cc_var_list, alpha_list]
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
n = length(x0);
xk = x0; % initializing
err = 10000; % err for convergence
k = 0; % counter
alphak = 100;  % initializing
t_x_list = zeros(n, kmax);
t_alpha_list = zeros(1, kmax);  % initializing convergence list
t_cc_grad_list = zeros(1, kmax);  % initializing convergence list
t_cc_func_list = zeros(1, kmax);  % initializing convergence list
t_cc_dgrad_list = zeros(1, kmax);  % initializing convergence list
t_cc_var_list = zeros(1, kmax);  % initializing convergence list
[grad_old, ~] = gradest(f,xk);  % initializing
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
    
    % RECORD
    for i = 1:n
        t_x_list(i, k+1) = xk(i);
    end
    t_alpha_list(k+1) = alphak; % record step
    t_cc_grad_list(k+1) = cc_grad(grad); % record convergence
    t_cc_func_list(k+1) = cc_func(f(xk_old), f(xk)); % record convergence
    t_cc_dgrad_list(k+1) = cc_dgrad(grad_old, grad); % record convergence
    t_cc_var_list(k+1) = cc_var(xk_old, xk); % record convergence
    
    % CONVERGENCE
    err = cc_grad(grad); % condition of gradient
    % err = cc_func(f(xk_old), f(xk)); % condition of objective function difference
    % err = cc_dgrad(grad_old, grad); % condition of gradient difference
    % err = cc_var(xk_old, xk)); % condition of parameters difference
    % err = min(cc_func(f(xk_old), f(xk)), min(cc_grad(grad), cc_var(xk_old, xk))); % choose value for convergence
    
    % COUNTUP
    k = k + 1; % count
    
end

%%
x_list = zeros(n, k + 1);  % initializing return x list
alpha_list = zeros(1, k);  % initializing return alpha list
cc_grad_list = zeros(1, k);  % initializing return convergence list
cc_func_list = zeros(1, k);  % initializing return convergence list
cc_dgrad_list = zeros(1, k);  % initializing return convergence list
cc_var_list = zeros(1, k);  % initializing return convergence list
for j = 1:n
    x_list(j, 1) = x0(j);
end
for i = 1:k
    for j = 1:n
        x_list(j, i+1) = t_x_list(j, i);
    end
    alpha_list(i) = t_alpha_list(i);
    cc_grad_list(i) = t_cc_grad_list(i);
    cc_func_list(i) = t_cc_func_list(i);
    cc_dgrad_list(i) = t_cc_dgrad_list(i);
    cc_var_list(i) = t_cc_var_list(i);
end
x = xk;

end


function main()
%f = @(x, p) (x(1)+3 - p(1))^6 + (x(2)+1 - p(2))^4;
%p = [1, 0.8548, 1, 1.3514];
p = [1, 2.11, 1, 2.11];
x0 = [2, 2];
f = @(x) energy_H2O(x, p);
x_obj = [4, 4];
p0 = [3, 3];
x00 = [1, 1];
tol = 0.000001;
kmax = 10000;
x_conj = conjugate_gradient(f, x0, tol, kmax);
x_stee = steepest_descent(f, x0, tol, kmax);
disp(x_conj)
disp(x_stee)
%[x_app, p_app, x_list, p_list, cc_grad_list, cc_func_list, cc_dgrad_list, cc_var_list, alpha_list] = conjugate_gradient_param_record(f, x00, p0, x_obj, tol, kmax);
%[x_app, p_app, x_list, p_list, cc_grad_list, cc_func_list, cc_dgrad_list, cc_var_list, alpha_list] = steepest_descent_param_record(f, x00, p0, x_obj, tol, kmax);
%disp(x_app); % 2.2449
%disp(p_app);
%{
% Probleme direct
p = [1, 0.8548, 1, 1.3514];
x0 = [4, 2];
f = @(x) energy_H2O(x, p);
df1 = @(x) da_energy_H2O(x, p);
df2 = @(x) dtheta_energy_H2O(x, p);
%df = {df1, df2};
%df1 = @(x) dljp(x, [1, 1]);
%df = {df1};
%[x_opt, x_list, cc_grad_list, cc_func_list, cc_dgrad_list, cc_var_list, alpha_list] = steepest_descent_record(f, x0, 0.0001, 100, df);
%[x_opt, x_list, cc_grad_list, cc_func_list, cc_dgrad_list, cc_var_list, alpha_list] = steepest_descent_record(f, x0, 0.0001, 100);
x_opt = conjugate_gradient(f, x0, 0.0001, 100);
disp(x_opt)
figure(1)
plot(log(cc_grad_list))
figure(2)
plot(log(cc_func_list))
figure(3)
plot(log(cc_dgrad_list))
figure(4)
plot(log(cc_var_list))
figure(5)
plot(log(alpha_list))
figure(6)
plot(x_list(1, :), x_list(2, :), 'bo-')
%}

% Probleme inverse
%{
f = @(x, p) energy_H2O(x, [1, p(1), 1, p(2)]);
x_obj = [0.9584, ]; % 
p0 = [3, 3];
x00 = [2, 1];
tol = 0.00001;
kmax = 1000;
[x_app, p_app, x_list, p_list, cc_grad_list, cc_func_list, cc_dgrad_list, cc_var_list, alpha_list] = steepest_descent_param_record_test(f, x00, p0, x_obj, tol, kmax);
disp(x_app); % 2.2449
disp(p_app);
figure(1)
plot(log(cc_grad_list))
figure(2)
plot(log(cc_func_list))
figure(3)
plot(log(cc_dgrad_list))
figure(4)
plot(log(cc_var_list))
figure(5)
plot(log(alpha_list))
figure(6)
plot(x_list(1, :), x_list(2, :), 'bo-')
figure(7)
plot(p_list(1, :), p_list(2, :), 'bo-')
%}
end

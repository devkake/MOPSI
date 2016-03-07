function test()
% Probleme inverse
f = @(x, p) energy_H2O(x, [1, p(1), 1, p(2)]);
x_obj = [0.9584, 1.8840]; % 
p0 = [2, 2];
x00 = [1.5, 1];
tol = 0.00001;
kmax = 1000;
[x_app, p_app, x_list, p_list, cc_grad_list, cc_func_list, cc_dgrad_list, cc_var_list, alpha_list] = steepest_descent_param_record_H2O(f, x00, p0, x_obj, tol, kmax);
%[x_app, p_app, x_list, p_list, cc_grad_list, cc_func_list, cc_dgrad_list, cc_var_list, alpha_list] = steepest_descent_param_record(f, x00, p0, x_obj, tol, kmax);
disp(x_app); % 2.2449
disp(p_app);
figure(1)
plot(log10(cc_grad_list))
figure(2)
plot(log10(cc_func_list))
figure(3)
plot(log10(cc_dgrad_list))
figure(4)
plot(log10(cc_var_list))
figure(5)
plot(log10(alpha_list))
figure(6)
plot(x_list(1, :), x_list(2, :), 'bo-')
figure(7)
plot(p_list(1, :), p_list(2, :), 'bo-')

end

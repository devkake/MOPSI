function result_CH4()
f = @(x, p) energy_CH4(x, [1, p(1), 1, p(2)]);
x_obj = [1.0870, 1.9021];  %  1.0870 A / 109.00 
p0 = [2, 2];
x00 = [2, 2];
tol = 0.0001;
kmax = 1000;
[x_app, p_app, x_list, p_list, cc_grad_list, cc_func_list, cc_dgrad_list, cc_var_list, alpha_list] = steepest_descent_param_record(f, x00, p0, x_obj, tol, kmax);
disp(x_app);
disp(p_app);
figure(1)
plot(log10(cc_grad_list))
title('Convergence of gradient')
xlabel('Iteration')
ylabel('Log10( convergence value )')
figure(2)
plot(log10(cc_func_list))
title('Convergence of difference of objective function')
xlabel('Iteration')
ylabel('Log10( convergence value )')
figure(3)
plot(log10(cc_dgrad_list))
title('Convergence of deifference of gradient')
xlabel('Iteration')
ylabel('Log10( convergence value )')
figure(4)
plot(log10(cc_var_list))
title('Convergence of deifference of variable')
xlabel('Iteration')
ylabel('Log10( convergence value )')
figure(5)
plot(log10(alpha_list))
title('Step')
xlabel('Iteration')
ylabel('Log10( step )')
figure(6)
plot(x_list(1, :), x_list(2, :), 'bo-')
title('Transition of configurational variables')
xlabel('$$a\ [\AA]$$', 'interpreter', 'latex','FontSize',16 )
ylabel('$$\theta\ [rad]$$', 'interpreter', 'latex','FontSize',16 )
figure(7)
plot(p_list(1, :), p_list(2, :), 'bo-')
title('Transition of parameters')
xlabel('$$\sigma_{CH}\ [\AA]$$', 'interpreter', 'latex','FontSize',16 )
ylabel('$$\sigma_{HH}\ [\AA]$$', 'interpreter', 'latex','FontSize',16 )
end
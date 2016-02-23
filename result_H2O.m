function result_H2O()
f = @(x, p) energy_H2O(x, [1, p(1), 1, p(2)]);
x_obj = [0.9584, 1.8230];  %  0.9584 A / 104.45 
p0 = [2, 2];
x00 = [2, 2];
tol = 0.0001;
kmax = 1000;
%[x_app, p_app, x_list, p_list, cc_grad_list, cc_func_list, cc_dgrad_list, cc_var_list, alpha_list] = steepest_descent_param_record(f, x00, p0, x_obj, tol, kmax);
[x_app, p_app, x_list, p_list, cc_grad_list, cc_func_list, cc_dgrad_list, cc_var_list, alpha_list] = conjugate_gradient_param_record(f, x00, p0, x_obj, tol, kmax);
disp(x_app); % 2.2449
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
xlabel('$$\sigma_{OH}\ [\AA]$$', 'interpreter', 'latex','FontSize',16 )
ylabel('$$\sigma_{HH}\ [\AA]$$', 'interpreter', 'latex','FontSize',16 )
%{
figure(8)
plot(p_list(1, :), p_list(2, :), 'bo-')
title('Transition of parameters')
xlabel('$$\sigma_{OH}\ [\AA]$$', 'interpreter', 'latex','FontSize',16 )
ylabel('$$\sigma_{HH}\ [\AA]$$', 'interpreter', 'latex','FontSize',16 )
hold on
sigma_oh = linspace(min(p_list(1, :)) - 0.1, max(p_list(1, :)) + 0.1, 10);
sigma_hh = linspace(min(p_list(2, :)) - 0.1, max(p_list(2, :)) + 0.1, 10);
[OH, HH] = meshgrid(sigma_oh, sigma_hh);
noh = length(sigma_oh);
nhh = length(sigma_hh);
Z = zeros(noh, nhh);
for i = 1:noh
    for j = 1:nhh
        Z(i,j) = objective_function(f, [sigma_oh(i), sigma_hh(j)], [2, 2], x_obj, tol, 1000);
    end
end
%f_dist = @(x,y) x^2 + y^2;
%f_dist = @(oh, hh) (2^(1/6)*oh - x_obj(1))^2 + (acos(1-((2^(1/6)*hh)^2*0.5/(2^(1/6)*oh)^2)) - x_obj(2))^2;
%disp(f_dist)
contour(OH, HH, Z, 10)
%}
end
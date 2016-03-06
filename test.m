function test()
x = [1, 1];
p = [1, 1, 1, 1];
n = length(x);
m = length(p);
disp(hessian_xx_energy_H2O(x, p));
disp(hessian_xp_energy_H2O(x, p));
fxp = @(x) energy_H2O(x(1:n), x(n+1:n+m));
xp = [x, p];
hess = hessian(fxp,xp);
Dxx = hess(1:n,1:n);
Dpx = hess(1:n,n+1:n+m);
grad_x_p = - Dxx \ Dpx;
disp(Dxx);
disp(Dpx);
disp(grad_x_p);
end

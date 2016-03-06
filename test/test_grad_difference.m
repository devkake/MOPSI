function test_grad_difference()
n = 50;
xl = linspace(0.1, 2.1, n);
yl = linspace(0.1, 2.1, n);
z = zeros(n, n);
p = [1, 1, 1, 1];
f = @(x) energy_H2O(x(1:2), x(3:6));
for i = 1:n
    for j = 1:n
        x = [xl(i), yl(j)];
        hxx_a = hessian_xx_energy_H2O(x, p);
        hxp_a = hessian_xp_energy_H2O(x, p);
        g_a = -hxx_a \ hxp_a;
        h_n = hessian(f,[x, p]);
        hxx_n = h_n(1:2,1:2);
        hxp_n = h_n(1:2,3:6);
        g_n = - hxx_n \ hxp_n;
        z(i, j) = log10(norm(g_a - g_n, 2)/norm(g_a, 2));
    end
end
figure(3);
surf(xl, yl, z);
end

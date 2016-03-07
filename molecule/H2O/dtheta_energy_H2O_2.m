function dftheta = dtheta_energy_H2O_2(x, p)
a = x(1);
theta = x(2);
b = a * sqrt(2 * (1 - cos(theta)));
p_hh = [1, p(2)];
dftheta = dljp(b, p_hh) * (sin(theta) / (2 * sqrt(1 - cos(theta))));
end


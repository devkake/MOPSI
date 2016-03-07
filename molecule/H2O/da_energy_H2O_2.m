function dfa = da_energy_H2O_2(x, p)
a = x(1);
theta = x(2);
b = a * sqrt(2 * (1 - cos(theta)));
p_oh = [1, p(1)];
p_hh = [1, p(2)];
dfa = 2*dljp(a, p_oh) + dljp(b, p_hh) * (sqrt(2 * (1 - cos(theta))));
end


function dfa = da_energy_H2O(x, p)
a = x(1);
theta = x(2);
b = a * sqrt(2 * (1 - cos(theta)));
p_oh = [p(1), p(2)];
p_hh = [p(3), p(4)];
dfa = 2*dljp(a, p_oh) + dljp(b, p_hh) * (sqrt(2 * (1 - cos(theta))));
end


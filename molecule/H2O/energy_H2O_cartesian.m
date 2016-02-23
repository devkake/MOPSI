function energy = energy_H2O_cartesian(x, p)
xh1 = x(1:3);
xo = x(4:6);
xh2 = x(7:9);
a1 = norm(xh1 - xo, 2);
a2 = norm(xh2 - xo, 2);
b = norm(xh1 - xh2, 2);
p_oh = [p(1), p(2)];
p_hh = [p(3), p(4)];
energy = ljp(a1, p_oh) + ljp(a2, p_oh) + ljp(b, p_hh);
end
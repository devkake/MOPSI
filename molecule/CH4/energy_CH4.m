function energy = energy_CH4(x, p)
% methane
% a = 1.0870 A
% theta = 109 / 180 * pi <- fixe
a = x;
theta = 1.9021;
b = a * sqrt(2 * (1 - cos(theta)));
p_ch = [p(1), p(2)];
p_hh = [p(3), p(4)];
energy = 4*ljp(a, p_ch) + 6*ljp(b, p_hh);

end
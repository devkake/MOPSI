function energy = energy_C2H6(x, p)
% ethane
% cc = 1.535 A
% ch = 1.092 A
% theta = 109 / 180 * pi (angle C-C-H)
cc = x(1);
ch1 = x(2);
theta = x(3);
ch2 = sqrt(cc^2 + ch1^2 - 2 * cc * ch1 * cos(theta));
hh1 = sqrt(3) * ch1 * cos(theta);
hh2 = sqrt(hh1^2 / 3 + (cc + 2 * ch1 * sin(theta))^2);
hh3 = sqrt(4 * hh1^2 / 3 + (cc + 2 * ch1 * sin(theta))^2);
p_cc = [p(1), p(2)];
p_ch = [p(3), p(4)];
p_hh = [p(5), p(6)];

energy = ljp(cc, p_cc) + 6*ljp(ch1, p_ch) + 6*ljp(ch2, p_ch) + 6*ljp(hh1, p_hh) + 6*ljp(hh2, p_hh) + 3*ljp(hh3, p_hh);

end
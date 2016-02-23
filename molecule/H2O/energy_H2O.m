function energy = energy_H2O(x, p)
% a = 0.9584 A
% theta = 104.45 / 180 * pi
a = x(1);
theta = x(2);
b = a * sqrt(2 * (1 - cos(theta)));
p_oh = [p(1), p(2)];
p_hh = [p(3), p(4)];
energy = 2*ljp(a, p_oh) + ljp(b, p_hh);

%{
a = conf(1);
theta = conf(2);
epsilon_OH = param(1); % Angstrom
epsilon_HH = param(2); % Angstrom
sigma_OH = param(3); % KJ / mol
sigma_HH = param(4); % KJ / mol
b = a * sqrt(2 * (1 - cos(theta)));
energy = 2 * lennard_jones_potential(a, epsilon_OH, sigma_OH) + lennard_jones_potential(b, epsilon_HH, sigma_HH);
a = x(1);
theta = x(2);
epsilon1 = 0.25; % KJ / mol
sigma1 = 2.930; % A
param1 = [epsilon1, sigma1];
epsilon2 = 0.651696; % KJ / mol
sigma2 = 3.16555; % A
param2 = [epsilon2, sigma2];
b = sqrt(2*(a^2)*(1-cos(theta)));
energy = ljp(b, param1) + 2*ljp(a, param2);
%}
end


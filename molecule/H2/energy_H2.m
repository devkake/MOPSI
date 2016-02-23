function energy = energy_H2(x, p)
% a = 0.7414 A
a = x(1);
p_hh = [p(1), p(2)];
energy = ljp(a, p_hh);
%{
epsilon = 37.0;  % kb K
sigma = 2.93; % A
param = [epsilon, sigma];
energy = ljp(a, param);
epsilon_HH = param(1); % Angstrom
sigma_HH = param(2); % KJ / mol
energy = lennard_jones_potential(a, epsilon_HH, sigma_HH);
%}
end

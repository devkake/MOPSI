function H = hessian_xp_energy_H2O(x, p)
b = sqrt(2)*x(1)*(1-cos(x(2)))^0.5;
p_oh = [1, p(1)];
p_hh = [1, p(2)];
H = zeros(2,2);
H(1, 1) = 2*d2_LJ_dx_d_sigma(x(1), p_oh);
H(1, 2) = d2_LJ_dx_d_sigma(b, p_hh) * d_b_d_a(x);
H(2, 1) = 0;
H(2, 2) = d2_LJ_dx_d_sigma(b, p_hh) * d_b_d_theta(x);
%{
b = sqrt(2)*x(1)*(1-cos(x(2)))^0.5;
p_oh = [p(1), p(2)];
p_hh = [p(3), p(4)];
H = zeros(2,4);
H(1, 1) = 2*d2_LJ_dx_d_epsilon(x(1), p_oh);
H(1, 2) = 2*d2_LJ_dx_d_sigma(x(1), p_oh);
H(1, 3) = d2_LJ_dx_d_epsilon(b, p_hh) * d_b_d_a(x);
H(1, 4) = d2_LJ_dx_d_sigma(b, p_hh) * d_b_d_a(x);
H(2, 1) = 0;
H(2, 2) = 0;
H(2, 3) = d2_LJ_dx_d_epsilon(b, p_hh) * d_b_d_theta(x);
H(2, 4) = d2_LJ_dx_d_sigma(b, p_hh) * d_b_d_theta(x);
%}
end

function y = d2_LJ_dx_d_epsilon(x, p)
y = 24*p(2)^6/x^7 - 48*p(2)^12/x^13;
end

function y = d2_LJ_dx_d_sigma(x, p)
y = 144*p(1)*p(2)^5/x^7 - 576*p(1)*p(2)^11/x^13;
end

function y = d_b_d_a(x)
y = sqrt(2)*(1-cos(x(2)))^0.5;
end

function y = d_b_d_theta(x)
y = sqrt(2)*x(1)*0.5*(1-cos(x(2)))^(-0.50)*(sin(x(2)));
end
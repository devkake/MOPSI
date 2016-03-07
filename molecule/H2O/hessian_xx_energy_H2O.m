function H = hessian_xx_energy_H2O(x, p)
b = sqrt(2)*x(1)*(1-cos(x(2)))^0.5;
% p_oh = [p(1), p(2)];
% p_hh = [p(3), p(4)];
p_oh = [1, p(1)];
p_hh = [1, p(2)];
H = zeros(2,2);
H(1, 1) = 2*d2_LJ_d_x2(x(1),p_oh) + d2_LJ_d_x2(b,p_hh)*(d_b_d_a(x))^2;
H(1, 2) = d2_LJ_d_x2(b,p_hh)*d_b_d_a(x)*d_b_d_theta(x) + d_LJ_d_x(b,p_hh)*(d2_b_d_a_dtheta(x));
H(2, 1) = H(1, 2);
H(2, 2) = d2_LJ_d_x2(b,p_hh)*d_b_d_theta(x)^2 + d_LJ_d_x(b,p_hh)*d2_b_d_theta2(x);
end

function y = d_LJ_d_x(x, p)
y = 24*p(1)*p(2)^6/x^7 - 48*p(1)*p(2)^12/x^13;
end

function y = d2_LJ_d_x2(x, p)
y = -168*p(1)*p(2)^6/x^8 + 624*p(1)*p(2)^12/x^14;
end

function y = d_b_d_a(x)
y = sqrt(2)*(1-cos(x(2)))^0.5;
end

function y = d_b_d_theta(x)
y = sqrt(2)*x(1)*0.5*(1-cos(x(2)))^(-0.50)*(sin(x(2)));
end

function y = d2_b_d_theta2(x)
y = sqrt(2)*x(1)*0.5*((-0.5*(1-cos(x(2)))^(-1.5)*(sin(x(2)))^2) + ((1-cos(x(2)))^(-0.50)*(cos(x(2)))));
end

function y = d2_b_d_a_dtheta(x)
y = sqrt(2)*0.5*(1-cos(x(2)))^(-0.5)*sin(x(2));
end
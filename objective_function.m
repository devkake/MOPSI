function d = objective_function(f, p, x0, x_obj, tol, kmax, df)
f_x = @(x) f(x, p);
% x_app = 0;
if (nargin == 7)
    dfa = @(x) df{1}(x, p);
    dftheta = @(x) df{2}(x, p);
    dfx = {dfa, dftheta};
    x_app = steepest_descent(f_x, x0, tol, kmax, dfx);
else
    x_app = steepest_descent(f_x, x0, tol, kmax);
    % x_app = conjugate_gradient(f_x, x0, tol, kmax); % error
end
d = norm(x_obj - x_app, 2)^2;
end
function [r] = cc_dgrad(gradient_before, gradient_after)
r = norm(gradient_after - gradient_before, 2) / norm(gradient_before, 2);
end

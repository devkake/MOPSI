function [r] = cc_var(variable_before, variable_after)
r = norm(variable_after - variable_before) / norm(variable_before);
end
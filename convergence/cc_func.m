function [r] = cc_func(function_before, function_after)
r = abs(function_after-function_before) / abs(function_before);
end
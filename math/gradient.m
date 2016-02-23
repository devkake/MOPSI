function [ grad ] =  gradient(f, x, h)
%GRADIENT Obtaine gradient of function with parameter h
%   completed? : FALSE
%   comment : name confliction
%   comment : h should be changed in proportion to x value
%   f : a function to differentiate
%   x : a point to differentiate
%   h : difference for differential

n = length(x);
grad = zeros(n,1);
for i = 1 : n
    grad(i) =  (f(x(i) + h) - f(x(i)))/h;
end
%grad = (f(x. + h) - f(x.)) / h;
end


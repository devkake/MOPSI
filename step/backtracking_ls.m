function [alpha] = backtracking_ls(f, x, p, grad_f, alpha0, tau, c, jmax)
%BACKTRACKING_LS backtracking line search
%   completed? : FALSE
%   alpha0 : initial alpha
%   f : objective function to minimize
%   x : point to differentiate
%   p : direction to search 
%   grad_f : gradient at df/dx at x
%   jmax : maximum iteration

%%
m = p'*grad_f; % direction
t = -c * m; % direction inclinated
%disp(t)
alphak = alpha0; % initializing
j = 0; % counter
y = f(x);
while(y - f(x + alphak*p') < alphak * t && j < jmax)
    j = j + 1; % counter
    alphak = alphak * tau; % 
end
alpha = alphak;

end
        
  
        
        
function [x] = conjugate_gradient(f, x0, tol, kmax)


err = 100; % err for convergence
k = 0; % counter
alpha0 = 1;  % initializing

[grad, ~] = gradest(f,x0);  % initializing
pk = -grad; % initializing
sk = pk;
alphak =alpha0;
xk = x0;
% alphak = backtracking_ls(alpha0, f, x0, pk, grad0, 1000);
% xk = x0 + alphak * pk;
% [grad, ~] = gradest(f,xk);

%%

while (k < kmax && err > tol)
    
    xk_old = xk; % record former position
    pk_old = pk;
    sk_old = sk;
    %grad_old = grad;
    
    
    [grad, ~] = gradest(f,xk);
    pk = -grad;
    betak = (pk*pk')/(pk_old*pk_old');
    %betak = (pk*(pk'-pk_old'))/(pk*pk');
    %betak = max(betak, 0);
    sk = pk + betak * sk_old;
    alphak = backtracking_ls(f, xk, sk', grad', alphak*2, 0.5, 0.5, kmax);
    xk = xk_old + alphak*sk; % update position

    err = cc_grad(grad);
    % err = cc_var(xk_old, xk);
    
    k = k + 1; % count
    
end

x = xk;

end
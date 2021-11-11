function [x,a, g, err, f, stopiter,X] = Newton(x,fun,gfun, Hfun,xtrue)
% Newton update 
% set params and allocate memory
tol = 1e-5;
iter =1 ;
itermax = 1e4;
err = zeros(itermax,1);
f = zeros(itermax,1);
g = zeros(itermax,1);
a = zeros(itermax,1);
X = zeros(itermax,2); 
alpha =1; 
% initial guess
f(1) = fun(x); 
g(1) = norm(gfun(x)); 
grad = gfun(x); 
err(1) = norm(x-xtrue); 
a(1) = alpha; 
X(1,:) = x; 
while norm(grad) >tol && err(iter)>tol && iter < itermax 
    p = -Hfun(x)\grad; 
    alpha = linesearch(x,p,grad,alpha, fun, f(iter));
    x(:) =x +alpha*p; 
    grad =gfun(x); 
    X(iter+1,:) = x; 
    g(iter+1) = norm(grad); 
    f(iter+1) = fun(x); 
    err(iter+1) = norm(x-xtrue); 
    a(iter+1) = alpha; 
    iter =iter+1; 
end 
stopiter = iter; 
fprintf("Stop at iteration %d\n", stopiter)


end


function [alpha] = linesearch(x,p,g,alpha,fun,f0)
rho = 0.9; 
c = 1e-4; 
while fun(x+alpha*p) > f0+c*alpha*g'*p 
    alpha = rho*alpha; 
end
end


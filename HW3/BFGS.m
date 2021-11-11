function [x,a, g, err, f, stopiter, X] = BFGS(x,fun,gfun,xtrue)
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
Hi = eye(2);
gradp = zeros(size(grad));
xp = zeros(size(x));
X(1,:) = x; 
while norm(grad) >tol && err(iter)>tol && iter < itermax 
    y = grad - gradp; 
    s = x-xp;
    rho = 1/(y'*s);
    Hi = (eye(2)-rho*s*y')*Hi*(eye(2)-rho*y*s')+rho*s*s'; % inverse hessian approximation
    p = -Hi*grad; 
    alpha = linesearch(x,p,grad,alpha, fun, f(iter));
    gradp = grad; 
    xp = x; 
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
rho = 0.2; 
c = 0.2; 
while fun(x+alpha*p) > f0+c*alpha*g'*p 
    alpha = rho*alpha; 
end
end


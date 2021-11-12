function [x,fa, normg,err, f, iter,X] = SteepestDescent(x,g,fun,xtrue)

tol = 1e-5; 
itermax =1e5; 
iter =1; 
normg = zeros(itermax,1); 
err  = zeros(size(normg));
f =  zeros(size(normg));
fa = zeros(size(err));
X = zeros(itermax,2);
grad = g(x);
normgrad =norm(grad); 
normg(1) = normgrad; 
normerr = norm(x-xtrue);
err(1) = normerr;
a = 1;
fa(1) = a; 
f(1) = fun(x);
% fprintf("iteration: %d, err: %d, normgrad: %d\n", iter, err(iter), normg(iter)); 
X(1,:) =x; 

while normgrad>tol && normerr>tol && iter <itermax
    
    a = linesearch(x,-grad,grad,a,fun,f(iter));
    fa(iter+1) = a; 
    x(:) = x- a*grad; 
    X(iter+1,:)=x;
    grad = g(x);
    normgrad =norm(grad); 
    normg(iter+1) = normgrad;
    normerr = norm(x-xtrue);
    err(iter+1) =normerr;
    f(iter+1) = fun(x);
    iter = iter+1 ; 
%     if mod(iter,1000)==0
%     fprintf("iteration: %d, err: %d, normgrad: %d\n", iter, err(iter), normg(iter));
%     end
end
stopiter = iter; 
fprintf("Stop at iteration %d\n", stopiter)
end

function [a] = linesearch(x,p,g,a,fun,f0)
rho = 0.1; 
c = 0.01; 
while fun(x+a*p) > f0+c*a*g'*p 
    a = rho*a; 
end
end

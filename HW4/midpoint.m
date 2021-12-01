function [ts,y] = midpoint(fun,y0,N,T)
% midpoint Euler for automous system 
% fun - f(y)
% y0 - initial condition
% N - number of steps 
% T - end time (assume starting time is 0)

h = T/N; 
y = zeros(N+1,4); 
ts = zeros(N+1,1); 
y(1,:) = y0; 
ts(1) = 0; 

for i=2:N+1
    m = y(i-1,:)+h/2*fun(y(i-1,:))';
    y(i,:) = y(i-1,:)+h*fun(m)';
    ts(i) = (i-1)*h;
    
end


end



function [ts,y] = AB3(fun,y0,N,T)
% 3rd order 3-step Adam-Bashforth method for automous ODE system 
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
   if i<4
    k1 = fun(y(i-1,:))';
    k2 = fun(y(i-1,:)+h/2*k1)';
    k3 = fun(y(i-1,:)+h/2*k2)';
    k4 = fun(y(i-1,:)+h*k3)';
    y(i,:) = y(i-1,:) +h*(1/6*k1+1/3*k2+1/3*k3+1/6*k4);
    ts(i) = (i-1)*h;
 
   else
        y(i,:) = y(i-1,:)+h*(23*fun(y(i-1,:))-16*fun(y(i-2,:))+5*fun(y(i-3,:)))'/12;
   end
   ts(i) = (i-1)*h;

end


end



fun = @(y)  [y(3) y(4)  -y(1)/(y(1)^2+y(2)^2)  -y(2)/(y(1)^2+y(2)^2) ]';
% fun1 = @(t,y)  [y(3) y(4)  -y(1)/(y(1)^2+y(2)^2)  -y(2)/(y(1)^2+y(2)^2) ]';
f0 = @(t) [cos(t) sin(t) -sin(t) cos(t)];
y0 = [1 0 0 1]';
methods = {'euler','trapozoidal','midpoint','Rk4','AB3'};
t1 = zeros(5,8); % pre-allocate memory to measure CPU time
e =  zeros(5,8);

%% problem 1-section b
T = 4*pi;
N = [100, 200, 300, 400, 500, 600, 700, 800];
h0 = T./N;
H = [h0; h0.^2; h0.^2; h0.^4; h0.^3];

tic, [ts, y] =euler(fun, y0, 100, T); toc

for i=1:length(methods)
    for j =1:length(N)
        m = char(methods(i));
        fh = str2func(m);
        y = zeros(N(j)+1,4);
        ts =zeros(N(j)+1,1);
        tic,
        [ts(:), y(:)] = fh(fun,y0,N(j),T);
        tmp = toc;
        t1(i,j)=tmp/N(j); 
        err = max(abs(f0(ts)-y),[],2);
        e(i,j)  = err(end);
    end
    
end
% estimate constant C for each method  
C= mean(e./H,2);
% estimate CPU time per step for each method 
tps = mean(t1,2);

p = [1 2 2 4 3];
tol = [1e-3, 1e-6, 1e-12]; 
Nt = zeros(5,3); 
for i=1:5
  for j =1:3
      Nt(i,j)= T*nthroot(C(i)/tol(j), p(i));
  end

end
tt = zeros(5,3); 
for i=1:3
   tt(:,i) = Nt(:,i).*tps;
end
%% section c
h = 0.01*pi; 
N = 3*1e3; % number of steps for AB3 to get error close to 0.1
% N = 7*1e5; % number of steps for RK4 to get error close to 0.1

% % test for RK4
% y = y0;
% for i=2:N+1
%     k1 = fun(y);
%     k2 = fun(y+h/2*k1);
%     k3 = fun(y+h/2*k2);
%     k4 = fun(y+h*k3);
%     y = y +h*(1/6*k1+1/3*k2+1/3*k3+1/6*k4);
%     ts = (i-1)*h;
%     
% end
% max(abs(f0(ts)'-y))
% test for AB3 
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
err = max(abs(f0(ts)-y),[],2);

err(end)


%% problem 2
% mu = 1e2; %1e2, 1e4, 1e6
tol = logspace(-9,-3,3);
y0 = [2; 0];
tspan = [0 200];
t= zeros(length(tol),3);

for i=1:length(tol)

options = odeset('Events',@myEventsFcn ,'AbsTol', tol(i),'RelTol', tol(i));
tic,
[t1,y1, te, ye, ie] = ode45(@func,tspan, y0, options); % ode45
t(i,1) = toc; 

tic, 
[t2,y2,te,ye,ie] = ode15s(@func,tspan,y0,options); % ode15s
t(i,2) = toc; 

tic, 
[t3,y3,te,ye,ie] = ode23s(@func,tspan,y0,options); % ode23s
t(i,3) = toc; 
end


figure()
loglog(tol, t(:,1)) % ode45
hold on
loglog(tol, t(:,2)) % ode 15s
hold on 
loglog(tol, t(:,3)) % ode 23s
hold off
legend("ode45","ode15s","ode23s")
xlabel("tolerance")
ylabel("CPU time")
title("log(CPU time) vs log(tol) for mu = 10")


function [value,isterminal,direction] = myEventsFcn(t,y)
value = y(1)-2;
isterminal = 1; 
direction = 1;

end
function dydt = func(t,y)
dydt = [y(2); 10*((1-y(1)^2)*y(2))-y(1)];

end

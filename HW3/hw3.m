%% Rosenbrock function 

fun = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2; 
gfun = @(x) [400*(x(1)^3-x(1)*x(2))-2*(1-x(1)) 200*(x(2)-x(1)^2)]';
Hfun = @(x) [1200*x(1)^2-400*x(2)+2 -400*x(1); -400*x(1) 200];
x0 = [1.2 1]'; 
xtrue = [1 1]'; 

[x,a1, normg,err1,f, iter1, XP] = SP(x0,gfun,fun,xtrue); % steepest descent 
[x,a2, g, err2, f, iter2, XN] = Newton(x0,fun,gfun, Hfun,xtrue);
[x,a3, g, err3, f, iter3, XS] = SR1(x0,fun,gfun,xtrue);
[x,a4, g, err4, f, iter4,XX] = BFGS(x0,fun,gfun,xtrue);

%% level curves
figure(1)
x = linspace(0.5,1.5);
y = linspace(0.5,1.5);
[X,Y] = meshgrid(x,y);
Z =  100*(Y-X.^2).^2+(1-X).^2; 
contour(X,Y,Z,30)
hold on 
plot(XX(1:iter4,1),XX(1:iter4,2), '--or') % 314 is the stopping iteration for BFGS
hold on 
plot(XS(1:iter3,1),XS(1:iter3,2), '-*b') % 144 is the stopping iteration for SR1
hold on 
plot(XN(1:iter2,1),XN(1:iter2,2), '--pk') % 107 is the stopping iteration for Newton
hold on 
plot(XP(1:iter1,1),XP(1:iter1,2), '-og') % 20084 is the stopping iteration for Steepest Descent
hold off
legend("","BFGS","SR1","Newton","SD") % SD stands for steepest descent
title("Level curve of Rosenbrock function with with x0=(1.2,1)")
xlabel("x values")
ylabel("y values")

%% plotting errors 
figure(2)
LS = 3;
semilogy(err1(1:iter1),'k','LineWidth',LS) %SD 
hold on 
% legend('SG')
% xlabel('iteration')
% ylabel('error norm')
% title('Error norm ||xtrue-x|| against iteration number')
% figure(3)
semilogy(err2(1:iter2),'b','LineWidth',LS) % Newton 
% legend('Newton')
% xlabel('iteration')
% ylabel('error norm')
% title('Error norm ||xtrue-x|| against iteration number')
hold on
semilogy(err3(1:iter3),'r','LineWidth',LS) % SR1
% legend('SR1')
% xlabel('iteration')
% ylabel('error norm')
% title('Error norm ||xtrue-x|| against iteration number')
hold on
semilogy(err4(1:iter4),'g','LineWidth',LS) % BFGS 
legend('SG','Newton','SR1','BFGS')
xlabel('iteration')
ylabel('error norm')
title('Error norm ||xtrue-x|| against iteration number with x0=(1.2,1)')

%% step size 
LS = 3;
figure(3)
% SD 
plot(a1(1:iter1),'LineWidth',LS)
hold on
% Newton 
plot(a2(1:iter2),'LineWidth',LS)
hold on
% SR1
plot(a3(1:iter3),'LineWidth',LS)
hold on
% BFGS 
plot(a4(1:iter4),'LineWidth',LS)
hold off
xlabel("iteration")
ylabel("step size")
title("Step size against iteration with x0=(1.2,1)")
legend("SD","Newton","SR1","BFGS")

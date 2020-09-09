close all; clear; clc
kx = linspace(0,pi/2,100);
Linear_Extra = (-1+2*exp(1i*kx)-exp(1i*(kx)))./(1+2*exp(-1i*kx)+exp(-2*1i*kx));
Quadra_Extra = (-1+3*exp(1i*kx)-3*exp(2*1i*kx)+exp(3*1i*kx))./...
    (1+3*exp(-1i*kx)+3*exp(-2*1i*kx)+exp(-3*1i*kx));
Zero_Grad = (-1+exp(1i*kx))./(1+exp(-1i*kx));
Antisym = (-1-exp(1i*kx))./(1-exp(-1i*kx));
Upwind_1st = (-1+exp(1i*kx)-1i*sin(kx))./(1+exp(-1i*kx)+1i*sin(kx));
Upwind_2nd = (-3+4*exp(1i*kx)-exp(2*1i*kx)-2*1i*sin(kx))./...
    (3+4*exp(-1i*kx)+exp(-2*1i*kx)+2*1i*sin(kx));
figure; hold on;
plot(kx,abs(Linear_Extra),'r--', kx,abs(Quadra_Extra),'k','LineWidth',2);
plot(kx,abs(Zero_Grad),'b',kx,abs(Antisym),'k*',kx,abs(Upwind_1st),'m-',kx,abs(Upwind_2nd),'go');
plot(kx,ones(size(kx)))
legend('Linear_Extrapolating','Quadratic Extrapolating','Homogenous Neumann',...
    'Antisymmetric','1st Order Upwind','2nd_Order Upwind','u_{N} = 0 ');
ylim([0 1.2])
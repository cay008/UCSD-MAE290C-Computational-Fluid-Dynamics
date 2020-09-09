close all;clear;clc;

%Well Done Padded Burgers
nu = 0.001; %Kinematic Viscosity;
CFL = 0.1; %CFL Number;
N = 128; %Original Grid
M = 192; %Extended Grid Points
dx = 2*pi/M; %Spatial Resolution;
x = 0:dx:(2*pi-dx); %Spatial Domain [0,2*pi);
x = x'; %Transfer it into Column Vector;


u = sin(x); %Initial Condition;
u_hat(:,1) = fft(u);%Fourier Coefficients of Initial Condition;

%Paramters of RK-Theta;
alpha1 = sqrt(2)-1; alpha2 = 0;  beta1 = 1-sqrt(2)/2; beta2 = beta1;
gamma1 = sqrt(2)/2; gamma2 = gamma1; zeta1 = 1-sqrt(2);

%Fourier Modes we are working with;
n = linspace(-M/2,M/2-1,M);
kn = reorder(n'); %Only For Function with Period of 2*pi, and change it into Column Vector;

figure;
i = 1; t = 0; dt = 0 ;
while t < 2
    %Adaptive Time Stepping
    dt = (CFL*dx)/abs((max(u(:,i))));
    
    %Zero Padding
    uv_hat = zeropadding(u_hat,N,M);
    
    %Application of RK-theta to advance fourier coefficient forward;
    u_hat = (u_hat - dt*alpha1*nu*(kn.^2).*u_hat - dt*gamma1*uv_hat)./...
        (1+dt*beta1*nu*kn.^2); %Here u_hat is essentially u_hat_star
    
    %Zero Padding
    u_star_hat = u_hat;
    uv_star_hat = zeropadding(u_star_hat,N,M);
    
    
    u_hat = (u_hat - dt*alpha2*nu*(kn.^2).*u_hat - dt*gamma2*uv_star_hat - dt*zeta1*uv_hat)./...
        (1 + dt*nu*beta2*kn.^2);
    
    %Inverse Fourier Transform to to back to physical domain;
    u(:,i+1) = ifft(u_hat);
    
    t = t + dt;
    plot(x,real(u(:,i+1))); title(num2str(t)); ylim([-1.5 1.5]);
    pause(0.0001);
    i = i + 1;
    
end


% Exchange the Order of Our Modes with that Given by Matlab fft
function u_fft = reorder(u_fft)
N = max(size(u_fft));
u_fft = -flip(u_fft);
%Store First N/2-1 elements
u_fft_inter = u_fft(1:(N/2-1));

u_fft(1:(N/2+1)) = u_fft((N/2):N);
u_fft((N/2+2):N) = u_fft_inter;
end
function uv_hat_0 = zeropadding(u_hat,N,M)
n = linspace(-M/2,M/2-1,M);
kn = reorder(n)';
u_hat_pad = u_hat ; u_hat_pad(N/2+2:M+1-N/2) = 0 ;
v_hat_pad = 1i*kn.*u_hat ; v_hat_pad(N/2+2:M+1-N/2) = 0;
uv_phy = ifft(u_hat_pad).*ifft(v_hat_pad);
uv_hat_0 = fft(uv_phy); uv_hat_0(N/2+2:M+1-N/2) = 0 ;
end


close all;clear;clc;
% Raw Burger
nu = 0.001; %Kinematic Viscosity;
CFL = 0.1; %CFL Number;
N = 128; %Total number of Grid Points;
dx = 2*pi/N; %Spatial Resolution;

x = 0:dx:(2*pi-dx); %Spatial Domain [0,2*pi);
x = x'; %Transfer it into Column Vector;

u = sin(x); %Initial Condition;
u_hat(:,1) = fft(u);%Fourier Coefficients of Initial Condition;

%Paramters of RK-Theta;
alpha1 = sqrt(2)-1; alpha2 = 0;  beta1 = 1-sqrt(2)/2; beta2 = beta1;
gamma1 = sqrt(2)/2; gamma2 = gamma1; zeta1 = 1-sqrt(2);

%Fourier Modes we are working with;
n = linspace(-N/2,N/2-1,N);
kn = reorder(n'); %Only For Function with Period of 2*pi, and change it into Column Vector;

figure;
i = 1; t = 0; dt = 100;

while t < 2
    %Adaptive Time Stepping
    dt = (CFL*dx)/abs((max(u(:,i))));
    
    %Apply Pseduospectral Method to deal with Non-Linear Terms;
    uv_hat = PseduoSpectral(u_hat,kn);
    
    %Application of RK-theta to advance fourier coefficient forward;
    u_hat = (u_hat - dt*alpha1*nu*(kn.^2).*u_hat - dt*gamma1*uv_hat)./...
        (1+dt*beta1*nu*kn.^2); %u_hat here is essentially u_hat_star;
    
    uv_star_hat = PseduoSpectral(u_hat,kn);
    
    u_hat = (u_hat - dt*alpha2*nu*(kn.^2).*u_hat - dt*gamma2*uv_star_hat - dt*zeta1*uv_hat)./...
        (1 + dt*nu*beta2*kn.^2);
    
    %Inverse Fourier Transform to back to physical domain;
    u(:,i+1) = ifft(u_hat);
    
    plot(x,real(u(:,i+1))); title(num2str(t)); ylim([-1.5,1.5])
    pause(0.0001);
    
    t = t + dt;
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

function uv_hat = PseduoSpectral(u_hat,kn)

v_hat = 1i*kn.*u_hat;
uv_phy = ifft(u_hat).*ifft(v_hat);
uv_hat = fft(uv_phy);
end

close all;clear;clc;
tic
% Chose Which Initial Condition to Use; 1 is purely random initial
% condition, 2 is random initial condition with energy spectrum;
initial_condition = 1;

% Physical Parameters
Re = 5e4;
CFL  = 0.5;

% Parameters of RK3
alpha1 = 29/96; alpha2 = -3/40; alpha3 = 1/6;
beta1 = 37/160; beta2 = 5/24; beta3 = 1/6;
gamma1 = 8/15; gamma2 = 5/12; gamma3 = 3/4;
zeta1 = -17/60; zeta2 = -5/12;

% Discretization
Nx = 512*3; Ny = Nx; n = Nx*Ny;
dx = 2*pi/Nx; dy = 2*pi/Ny;
x = 0:dx:2*pi-dx; y = 0:dy:2*pi-dy;

% Wave Numbers
kx = [0:Nx/2-1 -Nx/2:-1]; kx = repmat(kx.',1,Ny); kx = gpuArray(kx);
ky = [0:Ny/2-1 -Ny/2:-1]; ky = repmat(ky,Nx,1); ky = gpuArray(ky);
k = sqrt(kx.^2 + ky.^2 );

% Operators
L = (1i*kx).^2 + (1i*ky).^2;

%% Initial Conditions

if initial_condition == 1 % Pure Random, w/o satisfying energy spectrum
    
    w_p = randn(Nx,Ny); w_f = fft2(w_p); w_f = gpuArray(w_f);
    psi_f = -w_f./L;
    psi_f(1,1) = 0;
    psi_f = gpuArray(psi_f); % GPU Array
    u_f = 1i*psi_f.*ky; v_f = -1i*psi_f.*kx;
    u_p = real(ifft2(u_f)); v_p = real(ifft2(v_f));
    ens(1) = sum(w_p.^2,'all');
    E(1) = sum(0.5*(u_p.^2 + v_p.^2),'all');
    E_frac(1) = 1;
    ens_frac(1) = 1;
    
elseif initial_condition == 2 % Satisfy Energy Spectrum
    
    E_k = 4.*pi./3.*(2.*pi.*k).^2./(1+(2.*pi.*k).^2).^(5/3);
    w_f = abs(sqrt( (2*pi*k).^2.*E_k )); w_f = gpuArray(w_f);
    phase = rand(Nx,Ny);
    
    w_p = real(ifft2(w_f.*exp(1i*2*pi*phase)));
    
    w_f = fft2(w_p);
    psi_f = -w_f./L;
    psi_f(1,1) = 0;
    psi_f = gpuArray(psi_f); % GPU Array
    u_f = 1i*psi_f.*ky; v_f = -1i*psi_f.*kx;
    u_p = real(ifft2(u_f)); v_p = real(ifft2(v_f));
    ens(1) = sum(w_p.^2,'all');
    E(1) = sum(0.5*(u_p.^2 + v_p.^2),'all');
    E_frac(1) = 1;
    ens_frac(1) = 1;
end

%% PreLocation and Send All Variables to GPU
x = gpuArray(x); y = gpuArray(y);
Nx = gpuArray(Nx); Ny = gpuArray(Ny); n = gpuArray(n);
alpha1 = gpuArray(alpha1);alpha2 = gpuArray(alpha2);alpha3 = gpuArray(alpha3);
beta1 = gpuArray(beta1); beta2 = gpuArray(beta2); beta3 = gpuArray(beta3);
gamma1 = gpuArray(gamma1);gamma2 = gpuArray(gamma2);gamma3 = gpuArray(gamma3);
zeta1 = gpuArray(zeta1);zeta2 = gpuArray(zeta2);
Re = gpuArray(Re); CFL = gpuArray(CFL);
dx = gpuArray(dx); dy = gpuArray(dy);
L = gpuArray(L);
ens = gpuArray(ens);
E = gpuArray(E);
% E_frac = gpuArray(E_frac);
ens_frac = gpuArray(ens_frac);
H_f = zeros(Nx,Ny,'gpuArray');
H_star_f = zeros(Nx,Ny,'gpuArray');
H_star_star_f = zeros(Nx,Ny,'gpuArray');
w_p_plot = zeros(Nx,Nx,5,'gpuArray');
transient_time = zeros(1,9,'gpuArray');
%%
% Some Parameters
T = 5000;  T = gpuArray(T);
dt = 0.001; i = 1; t(1) = 0;  t = gpuArray(t);

% Video
video = VideoWriter('2D turbulence 1536 gpu.avi');


% Plot of Initial Condition
figure(1); pcolor(w_p); colormap(jet); shading interp; colorbar
string = sprintf('Time Step %f',i); title(string);
% frame = getframe(gcf);
% writeVideo(video,frame);
open(video);
for i = 1: 30000
    % Adaptive Time Stepping
    U = sqrt(u_p.^2 + v_p.^2);
    u_max = max(max(abs(U)));
    dt = CFL*dx/(pi*u_max);
    
%     if dt >= 0.001
%         dt = 0.001
%     end

    %------------Apply RK3 to March Vorticity Equation Forward------------%
    
    %% Zero Padding
    
    dwdx_f = 1i*w_f.*kx; dwdy_f = 1i*w_f.*ky ;
    % Dealiase Non-Linear Terms in Vorticity Equation By Zero Padding
    u_f(Nx/3+1:2/3*Nx , :) = 0 ; u_f(:,Ny/3+1:2/3*Ny) = 0 ;
    v_f(Nx/3+1:2/3*Nx , :) = 0 ; v_f(:,Ny/3+1:2/3*Ny) = 0 ;
    dwdx_f(Nx/3+1:2/3*Nx , :) = 0; dwdx_f(:,Ny/3+1:2/3*Ny) = 0 ;
    dwdy_f(Nx/3+1:2/3*Nx , :) = 0; dwdy_f(: , Ny/3+1:2/3*Ny) = 0;
    % Back to Physical Space
    u_p = real(ifft2(u_f)); v_p = real(ifft2(v_f));
    dwdx_p = real(ifft2(dwdx_f)); dwdy_p = real(ifft2(dwdy_f));
    % Non Linear Term in Physical Domain
    H_p = u_p.*dwdx_p + v_p.*dwdy_p;
    % Non Linear Term in Fourier/Frequency Domain
    H_f = fft2(H_p);
    % Zero Padding
    H_f(Nx/3+1:2/3*Nx , :) = 0;
    H_f(: , Ny/3+1:2/3*Ny) = 0;
    H_f = -H_f;
    %% Step 1 of RKW3
    
    % Solve Vorticity Equation
    w_f = (w_f + 1/Re*L*dt*alpha1.*w_f + gamma1*dt*H_f)./(1-1/Re*L*dt*beta1);
    % Solve for Streamfunction
    psi_f = -w_f./L;
    psi_f(1,1) = 0;
    
    %% Zero Padding
    
    u_f = 1i*psi_f.*ky;
    v_f = -1i*psi_f.*kx;
    dwdx_f = 1i*w_f.*kx; dwdy_f = 1i*w_f.*ky ;
    u_f(Nx/3+1:2/3*Nx , :) = 0 ; u_f(:,Ny/3+1:2/3*Ny) = 0 ;
    v_f(Nx/3+1:2/3*Nx , :) = 0 ; v_f(:,Ny/3+1:2/3*Ny) = 0 ;
    dwdx_f(Nx/3+1:2/3*Nx , :) = 0; dwdx_f(:,Ny/3+1:2/3*Ny) = 0 ;
    dwdy_f(Nx/3+1:2/3*Nx , :) = 0; dwdy_f(: , Ny/3+1:2/3*Ny) = 0;
    u_p = real(ifft2(u_f)); v_p = real(ifft2(v_f));
    dwdx_p = real(ifft2(dwdx_f)); dwdy_p = real(ifft2(dwdy_f));
    H_star_p = u_p.*dwdx_p + v_p.*dwdy_p;
    H_star_f = fft2(H_star_p);
    H_star_f(Nx/3+1:2/3*Nx , :) = 0;
    H_star_f(: , Ny/3+1:2/3*Ny) = 0;
    H_star_f = -H_star_f;
    %% Step 2 of RKW3
    w_f = (w_f + 1/Re*L*dt*alpha2.*w_f + gamma2*dt*H_star_f + zeta1*dt*H_f)./ ...
        (1 - 1/Re*L*dt*beta2); % w_star_f here is essential w_star_star_f
    psi_f = -w_f./L;
    psi_f(1,1) = 0;
    
    %% Zero Padding
    
    u_f = 1i*psi_f.*ky;
    v_f = -1i*psi_f.*kx;
    dwdx_f = 1i*w_f.*kx; dwdy_f = 1i*w_f.*ky ;
    u_f(Nx/3+1:2/3*Nx , :) = 0 ; u_f(:,Ny/3+1:2/3*Ny) = 0 ;
    v_f(Nx/3+1:2/3*Nx , :) = 0 ; v_f(:,Ny/3+1:2/3*Ny) = 0 ;
    dwdx_f(Nx/3+1:2/3*Nx , :) = 0; dwdx_f(:,Ny/3+1:2/3*Ny) = 0 ;
    dwdy_f(Nx/3+1:2/3*Nx , :) = 0; dwdy_f(: , Ny/3+1:2/3*Ny) = 0;
    u_p = real(ifft2(u_f)); v_p = real(ifft2(v_f));
    dwdx_p = real(ifft2(dwdx_f)); dwdy_p = real(ifft2(dwdy_f));
    H_star_star_p = u_p.*dwdx_p + v_p.*dwdy_p;
    H_star_star_f = fft2(H_star_star_p);
    H_star_star_f(Nx/3+1:2/3*Nx , :) = 0;
    H_star_star_f(: , Ny/3+1:2/3*Ny) = 0;
    H_star_star_f = -H_star_star_f;
    %% Step 3 of RKW3
    w_f = (w_f + 1/Re*L*dt*alpha3.*w_f + gamma3*dt*H_star_star_f + zeta2*dt*H_star_f)./ ...
        (1 - 1/Re*L*dt*beta3);
    psi_f = -w_f./L;
    psi_f(1,1) = 0;
    
    w_p = real(ifft2(w_f));
    
    i = i + 1;
    t(i) = t(i-1) + dt;
    
    u_f = 1i*psi_f.*ky; v_f = -1i*psi_f.*kx;
    u_p = real(ifft2(u_f)); v_p = real(ifft2(v_f));
    
    %% Enstrophy && Energy
    ens(i) = sum(w_p.^2,'all');
    ens_frac = ens(i)./ens(1);
    E(i) = sum(0.5*(u_p.^2 + v_p.^2),'all');
    E_frac(i) = gather(E(i)./E(1));
    
    %% Plot (Vdieo Gen)
    if mod(i,50) == 0
        pcolor(w_p); colormap(jet); shading interp; colorbar
        drawnow
        string = sprintf('Energy Fraction %f',E_frac(i));
        title(string)
        frame = getframe(gcf);
        writeVideo(video,frame);
    end
    
    if abs(E_frac(i) - 0.9) <= 0.02
        w_p_plot(:,:,1) = w_p;
        transient_time(1) = t(i);
    elseif abs(E_frac(i) - 0.8) <= 0.02
        w_p_plot(:,:,2) = w_p;
        transient_time(2) = t(i);
    elseif abs(E_frac(i) - 0.7) <= 0.005
        w_p_plot(:,:,3) = w_p;
        transient_time(3) = t(i);
    elseif abs(E_frac(i) - 0.6) <= 0.005
        w_p_plot(:,:,4) = w_p;
        transient_time(4) = t(i);
    elseif abs(E_frac(i) - 0.5) <= 0.005
        w_p_plot(:,:,5) = w_p;
        transient_time(5) = t(i);
    elseif abs(E_frac(i) - 0.4) <= 0.005
        w_p_plot(:,:,6) = w_p;
        transient_time(6) = t(i);
    elseif abs(E_frac(i) - 0.3) <= 0.005
        w_p_plot(:,:,7) = w_p;
        transient_time(7) = t(i);
    elseif abs(E_frac(i) - 0.2) <= 0.02
        w_p_plot(:,:,8) = w_p;
        transient_time(8) = t(i);
    elseif abs(E_frac(i) - 0.1) <= 0.02
        w_p_plot(:,:,9) = w_p;
        transient_time(9) = t(i);
    end
    
    if E_frac(i-1) <= 0.1
        break
    end
    
end

close(video);
figure; semilogy(t,ens); xlabel('time step'); ylabel('Enstrophy');
title('Time Evolution of Enstrophy'); saveas(gcf,'enstrophy.bmp');
figure; semilogy(t,E); xlabel('time step'); ylabel('Energy');
title('Time Evolution of Total Energy'); saveas(gcf,'energy.bmp');

figure; myplot(w_p_plot(:,:,1)); title('2D Vorticity with Energy of 90%');
saveas(gcf,'90.bmp'); saveas(gcf,'90.pdf')
figure; myplot(w_p_plot(:,:,2)); title('2D Vorticity with Energy of 80%');
saveas(gcf,'80.bmp'); saveas(gcf,'80.pdf')
figure; myplot(w_p_plot(:,:,3)); title('2D Vorticity with Energy of 70%');
saveas(gcf,'70.bmp'); saveas(gcf,'70.pdf')
figure; myplot(w_p_plot(:,:,4)); title('2D Vorticity with Energy of 60%');
saveas(gcf,'60.bmp'); saveas(gcf,'60.pdf')
figure; myplot(w_p_plot(:,:,5)); title('2D Vorticity with Energy of 50%');
saveas(gcf,'50.bmp'); saveas(gcf,'50.pdf')
figure; myplot(w_p_plot(:,:,6)); title('2D Vorticity with Energy of 40%');
saveas(gcf,'40.bmp'); saveas(gcf,'40.pdf')
figure; myplot(w_p_plot(:,:,7)); title('2D Vorticity with Energy of 30%');
saveas(gcf,'30.bmp'); saveas(gcf,'30.pdf')
figure; myplot(w_p_plot(:,:,8)); title('2D Vorticity with Energy of 20%');
saveas(gcf,'20.bmp'); saveas(gcf,'20.pdf')
figure; myplot(w_p_plot(:,:,9)); title('2D Vorticity with Energy of 10%');
saveas(gcf,'10.bmp') ;saveas(gcf,'10.pdf')

toc

figure; plot(t,ens); xlabel('time step'); ylabel('Enstrophy');
title('Time Evolution of Enstrophy');
figure; semilogy(t,ens); xlabel('time step'); ylabel('log(Enstrophy)');
title('Time Evolution of Enstrophy');
figure; plot(t,E); xlabel('time step'); ylabel('Energy');
title('Time Evolution of Total Energy');
figure; semilogy(t,E); xlabel('time step'); ylabel('log(Energy)');
title('Time Evolution of Total Energy');


function myplot(w_p)
pcolor(w_p); colormap(jet); shading interp; colorbar
end
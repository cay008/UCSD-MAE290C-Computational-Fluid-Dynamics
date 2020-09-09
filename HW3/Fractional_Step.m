close all;clear;clc; method = 2;
run('setup');

%Spatial Discretizaiton
dx = 2*pi/Nx; dy = 2*pi/Ny;
x = 0 : dx : 2*pi-dx; y = 0 : dy : 2*pi-dy;

%Inital Condition
u = sin(x).'*sin(y); v = cos(x).'*cos(y);
u_hat = fft2(u); v_hat = fft2(v); 
u_hat = reshape(u_hat,[],1); v_hat = reshape(v_hat,[],1); %Assemble Columns into Single Column Vector;
U_hat = [ u_hat; v_hat]; %Size: 2n by 1 ;

%                      Fractional Step Method                       %

%Step 1 - Find Intermediate Solution; Solve A*U_hat_star = U_hat - H_hat*dt

%Non-Linear Term;
H_u = repmat(sin(2*x).'/2,1,Ny); H_v = repmat(-sin(2*y)/2,Nx,1);
H_u_hat = fft2(H_u); H_v_hat = fft2(H_v); %Nx by Ny; Nx by Ny;
H_hat = [reshape(H_u_hat,[],1);reshape(H_v_hat,[],1)]; %2n by 1

A = A(1:2*n,1:2*n); B = inv(A);
RHS = U_hat - H_hat*dt; 
U_hat_star = A\RHS;

%Back to Grid Points - Fourier
U_hat_star = reshape(U_hat_star,Nx,2*Ny); %Size: 32 by 64;

%Back to Physical Domain
u_star = ifft2(U_hat_star(1:Nx,1:Ny));      %Intermediate x velocity ;
v_star = ifft2(U_hat_star(1:Nx,Ny+1:2*Ny)); %Intermediate y velocity ;

%Plots
figure; 
subplot(1,2,1); contourf(u_star); title('Intermediate u Velocity'); axis square;
subplot(1,2,2); contourf(v_star); title('Intermediate v velocity'); axis square;

%Step2 - Pressure ; Solve (DBG)*P = D*U_hat_star/dt

D = D(2*n+1:3*n,1:2*n);
G = G(1:2*n,2*n+1:3*n);
DBG = D*B*G; %Size n by n ; 
DBG(1,1) = 1;%Remove Singularity

U_hat_star = reshape(U_hat_star,[],1); %Back to 2n by 1
RHS = D*U_hat_star/dt;
P_hat = DBG\RHS;

%Back to Grid Poins - Fourier
P_hat = reshape(P_hat,Nx,Ny);
%Back to Physical
P = ifft2(P_hat);
figure; contourf(real(P)); title('P_n+1');

%Step3 - Correction Step; Solve U_hat_n+1 = U_hat_star + (BG)P*dt

BG = B * G;
P_hat = reshape(P_hat,[],1);
U_hat_star = U_hat_star - BG*P_hat*dt;
%Match with Grid Point;
U_hat_star = reshape(U_hat_star,Nx,2*Ny);
%Back to Physical Domain;
u = ifft2(U_hat_star(1:Nx,1:Ny));
v = ifft2(U_hat_star(1:Nx,Ny+1:2*Ny));

%Plots
figure;
subplot(1,2,1); contourf(real(u)); title('u_n+1'); axis square;
subplot(1,2,2); contourf(real(v)); title('v_n+1'); axis square;















close all; clear;clc; method = 2;
run('setup');

%Spatial Discretizaiton
dx = 2*pi/Nx; dy = 2*pi/Ny;
x = 0 : dx : 2*pi-dx ; y = 0 : dy : 2*pi-dy;

%Inital Condition
u = sin(x).'*sin(y); v = cos(x).'*cos(y);
u_hat = fft2(u); v_hat = fft2(v); 
%Aseemble Columns
u_hat = reshape(u_hat,[],1); v_hat = reshape(v_hat,[],1); 
U_hat = [ u_hat; v_hat]; %Size: 2n by 1 ;

%                       Artificial Compressbility                  %

%Non-Linear Term;
H_u = repmat(sin(2*x).'/2,1,Ny); H_v = repmat(-sin(2*y)/2,Nx,1);
H_u_hat = fft2(H_u); H_v_hat = fft2(H_v); %Nx by Ny; Nx by Ny;
H_hat = [reshape(H_u_hat,[],1);reshape(H_v_hat,[],1)]; %2n by 1 ;

%Solve Compressible NS
RHS = [ U_hat - H_hat*dt ; zeros(n,1) ];
Solution = M_ac\RHS;
u_hat = Solution(1:n); v_hat = Solution(n+1:2*n); P_hat = Solution(2*n+1:3*n);

%Back to Grid - Fourier
u_hat = reshape(u_hat,Nx,Ny); v_hat = reshape(v_hat,Nx,Ny);P = reshape(P_hat,Nx,Ny);

%Back to Physical
u = ifft2(u_hat); v = ifft2(v_hat); P = ifft2(P);

%Plot
figure;
subplot(1,2,1); contourf(abs(u)); title('u_n+1'); axis square;
subplot(1,2,2); contourf(abs(v)); title('v_n+1'); axis square;
figure;
contourf(abs(P)); title('P_n+1');















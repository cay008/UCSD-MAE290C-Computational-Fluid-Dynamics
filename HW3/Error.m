close all;clear;clc;
lambda = linspace(10e-4,0.1,10); s1 = size(lambda,2);
dt = linspace(10e-6,0.1,10);       s2 = size(dt,2);

for ii = 1 : s1
    for jj = 1 :s2
        
        [u_fs,v_fs,P_fs] = FS(lambda(ii),dt(jj));
        [u_ac,v_ac,P_ac] = AC(lambda(ii),dt(jj));
        u_err = u_fs - u_ac; v_err = v_fs - v_ac; P_err = P_fs - P_ac;
        
        u_err_norm(ii,jj) = norm(u_err);
        v_err_norm(ii,jj) = norm(v_err);
        P_err_norm(ii,jj) = norm(P_err);
        
    end
end
% x axis: Lambda; y axis: u/v/P;
figure; plot(lambda, u_err_norm); title('Error of u');
legend(sprintf('t=%f',dt(1)),sprintf('t=%f',dt(2)),sprintf('t=%f',dt(3)),...
    sprintf('t=%f',dt(4)),sprintf('t=%f',dt(5)),sprintf('t=%f',dt(6)),...
    sprintf('t=%f',dt(7)),sprintf('t=%f',dt(8)),sprintf('t=%f',dt(9)),...
    sprintf('t=%f',dt(10)));
xlabel('\lambda')
figure; plot(lambda,v_err_norm); title('Error of v');
xlabel('\lambda')
figure; plot(lambda,P_err_norm); title('Error of P');
xlabel('\lambda')


function [u,v,P] = FS(lambda,dt)
%Set up
Nx = 32; Ny = 32; n = Nx*Ny; Re = 1000;
alpha = [ 0:Nx/2-1 -Nx/2:-1 ];
beta  = [ 0:Ny/2-1 -Ny/2:-1 ];
iy = repmat(1i*alpha',[1,Ny]);
iiy = repmat(-alpha'.^2,[1,Ny]);
G = sparse(1:n,2*n+1:3*n,repmat(1i*alpha,[1,Ny]),3*n,3*n) + ...
    sparse(n+1:2*n,2*n+1:3*n,reshape(iy.',[],1),3*n,3*n);
D = transpose(G);
I = diag(ones(2*n,2*n));
I = sparse(1:2*n,1:2*n,I,3*n,3*n);
L = D*G;
L(1:n,1:n) = L(2*n+1:3*n,2*n+1:3*n);
L(n+1:2*n,n+1:2*n) = L(2*n+1:3*n,2*n+1:3*n);
L(2*n+1:3*n,2*n+1:3*n) = 0;
A = I - (dt/Re)*L;
M = A + dt*G + D;
B = inv(A(1:2*n,1:2*n));
%Spatial Discretizaiton
dx = 2*pi/Nx; dy = 2*pi/Ny;
x = 0 : dx : 2*pi-dx; y = 0 : dy : 2*pi-dy;
%Inital Condition
u = sin(x).'*sin(y); v = cos(x).'*cos(y);
u_hat = fft2(u); v_hat = fft2(v);
u_hat = reshape(u_hat,[],1); v_hat = reshape(v_hat,[],1); %Assemble Columns into Single Column Vector;
U_hat = [ u_hat; v_hat]; %Size: 2n by 1 ;

%Non-Linear Term;
H_u = repmat(sin(2*x).'/2,1,Ny); H_v = repmat(-sin(2*y)/2,Nx,1);
H_u_hat = fft2(H_u); H_v_hat = fft2(H_v); %Nx by Ny; Nx by Ny;
H_hat = [reshape(H_u_hat,[],1);reshape(H_v_hat,[],1)]; % n by n ;

%Fractional Step Method

%Step 1 - Find Intermediate Solution; Solve A*U_hat_star = U_hat - H_hat*dt
A = A(1:2*n,1:2*n);
RHS = U_hat - H_hat*dt;
U_hat_star = A\RHS;
%Back to Grid Points
U_hat_star = reshape(U_hat_star,Nx,2*Ny); %Size: 32 by 64;
%Back to Physical Domain
u_star = ifft2(U_hat_star(1:Nx,1:Ny));      %x velocity
v_star = ifft2(U_hat_star(1:Nx,Ny+1:2*Ny)); %y velocity

%Step2 - Pressure ; Solve (DBG)*P = D*U_hat_star/dt
D = D(2*n+1:3*n,1:2*n);
G = G(1:2*n,2*n+1:3*n);
DBG = D*B*G; %Size n by n ;
DBG(1,1) = 1;%Remove Singularity

U_hat_star = reshape(U_hat_star,[],1); %Back to 2n by 1
RHS = D*U_hat_star/dt;
P_hat = DBG\RHS;

%Back to Grid Poins
P_hat = reshape(P_hat,Nx,Ny);
%Back to Physical
P = ifft2(P_hat);

%Step3 - Correction Step; Solve U_hat_n+1 = U_hat_star + (BG)P*dt
BG = B * G;
P_hat = reshape(P_hat,[],1);
U_hat_star = U_hat_star + BG*P_hat*dt;
%Match with Grid Point;
U_hat_star = reshape(U_hat_star,Nx,2*Ny);
%Back to Physical Domain;
u = ifft2(U_hat_star(1:Nx,1:Ny));
v = ifft2(U_hat_star(1:Nx,Ny+1:2*Ny));
end
function [u,v,P] = AC(lambda,dt)
% Setup
Nx = 32; Ny = 32; n = Nx*Ny; Re = 1000;
alpha = [ 0:Nx/2-1 -Nx/2:-1 ];
beta  = [ 0:Ny/2-1 -Ny/2:-1 ];
iy = repmat(1i*alpha',[1,Ny]);
iiy = repmat(-alpha'.^2,[1,Ny]);
G = sparse(1:n,2*n+1:3*n,repmat(1i*alpha,[1,Ny]),3*n,3*n) + ...
    sparse(n+1:2*n,2*n+1:3*n,reshape(iy.',[],1),3*n,3*n);
D = transpose(G);
I = diag(ones(2*n,2*n));
I = sparse(1:2*n,1:2*n,I,3*n,3*n);
L = D*G;
L(1:n,1:n) = L(2*n+1:3*n,2*n+1:3*n);
L(n+1:2*n,n+1:2*n) = L(2*n+1:3*n,2*n+1:3*n);
L(2*n+1:3*n,2*n+1:3*n) = 0;
A = I - (dt/Re)*L;
M = A + dt*G + D;
I = diag(ones(n,n));
I = lambda*sparse(2*n+1:3*n,2*n+1:3*n,I,3*n,3*n);
M = M + I;
% Artificial Compressibility;
B = inv(A(1:2*n,1:2*n));
%Spatial Discretizaiton
dx = 2*pi/Nx; dy = 2*pi/Ny;
x = 0 : dx : 2*pi-dx; y = 0 : dy : 2*pi-dy;

%Inital Condition
u = sin(x).'*sin(y); v = cos(x).'*cos(y);
u_hat = fft2(u); v_hat = fft2(v);
u_hat = reshape(u_hat,[],1); v_hat = reshape(v_hat,[],1); %Assemble Columns into Single Column Vector;
U_hat = [ u_hat; v_hat]; %Size: 2n by 1 ;

%Non-Linear Term;
H_u = repmat(sin(2*x).'/2,1,Ny); H_v = repmat(-sin(2*y)/2,Nx,1);
H_u_hat = fft2(H_u); H_v_hat = fft2(H_v); %Nx by Ny; Nx by Ny;
H_hat = [reshape(H_u_hat,[],1);reshape(H_v_hat,[],1)]; % n by n ;

%Solve the Artifical Compressible NS
RHS = [U_hat-H_hat*dt;zeros(n,1)];
Solution = M\RHS;
u_hat = Solution(1:n); v_hat = Solution(n+1:2*n); P_hat = Solution(2*n+1:3*n);
%Match with Grid
u_hat = reshape(u_hat,Nx,Ny); v_hat = reshape(v_hat,Nx,Ny);P = reshape(P_hat,Nx,Ny);
%Back to Physical
u = ifft2(u_hat); v = ifft2(v_hat); P = ifft2(P);
end
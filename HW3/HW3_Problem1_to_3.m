close all;clear;clc; method = 2;
%Problem1
Nx = 32; Ny = 32; n = Nx*Ny; Re = 1000; dt = 0.01;
alpha = [ 0:Nx/2-1 -Nx/2:-1 ]; %Wave # in x direction
beta  = [ 0:Ny/2-1 -Ny/2:-1 ]; %Wave # in y direction

iy = repmat(1i*alpha',[1,Ny]);
iiy = repmat(-alpha'.^2,[1,Ny]);

G = sparse(1:n,2*n+1:3*n,repmat(1i*alpha,[1,Ny]),3*n,3*n) + ... 
    sparse(n+1:2*n,2*n+1:3*n,reshape(iy.',[],1),3*n,3*n);

D = transpose(G);

I = diag(ones(2*n,2*n));
I = sparse(1:2*n,1:2*n,I,3*n,3*n);

if method == 1;
    L = sparse(1:n,1:n,repmat(-alpha.^2,[1,Ny]),3*n,3*n) + ...
        sparse(1:n,1:n,reshape(iiy',1,[]),3*n,3*n) + ...
        sparse(n+1:2*n,n+1:2*n,repmat(-alpha.^2,[1,Ny]),3*n,3*n) + ...
        sparse(n+1:2*n,n+1:2*n,reshape(iiy',1,[]),3*n,3*n);
else
    L = D*G;
    L(1:n,1:n) = L(2*n+1:3*n,2*n+1:3*n);
    L(n+1:2*n,n+1:2*n) = L(2*n+1:3*n,2*n+1:3*n);
    %Back to 3n by 3n
    L(2*n+1:3*n,2*n+1:3*n) = 0;
end


A = I - (dt/Re)*L;

M = A + dt*G + D;

figure; imagesc(abs(M)); title('M')
figure; imagesc(abs(D)); title('D');
figure; imagesc(abs(G)); title('G');
figure; imagesc(abs(A)); title('A');


%Problem2
for i = 1:3*n
    if any(M(i,:)) == 0
        index = i %Index = 2049; i.e. (2n + 1) th row is full of 0s';
    end
end

%Problem3 
lambda = 0.01;
I = diag(ones(n,n));
I = lambda*sparse(2*n+1:3*n,2*n+1:3*n,I,3*n,3*n);

M_ac = M + I;
figure; imagesc(abs(M_ac)); title('Problem3 - New M')
ew = eigs(M_ac,1,'smallestabs')













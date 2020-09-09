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

L = D*G;
L(1:n,1:n) = L(2*n+1:3*n,2*n+1:3*n);
L(n+1:2*n,n+1:2*n) = L(2*n+1:3*n,2*n+1:3*n);
%Back to 3n by 3n
L(2*n+1:3*n,2*n+1:3*n) = 0;


A = I - (dt/Re)*L;

M = A + dt*G + D;

lambda = 0.01;
I = diag(ones(n,n));
I = lambda*sparse(2*n+1:3*n,2*n+1:3*n,I,3*n,3*n);

M_ac = M + I;
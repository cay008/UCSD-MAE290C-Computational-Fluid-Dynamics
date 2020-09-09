close all;clear;clc
% Option1: Run all the pronblems successivelly; [ 0 1 ]
% Option2: Run Specfic Problem; [ 1 2 3 4 5 6 7 ]
option = 1;
problem = 0;
% Setup
N = 200; dx = 10/N;
dt = 0.01;  T = 2; nof_t = T/dt;
x = 0:dx:10-dx;
c = (3/4)*sqrt(3)*dx/dt;
A = 3;

% Parameters of RKW3
a21 = 8/15; a32 = 5/12; b1 = 1/4; b3 = 3/4; b2 = 0;
zeta2 = b1 - a21; zeta3 = b2 - a32;

% Define Discretized Version of RHS convection term -c(du/dx)
RHS = @(right, left) (-c/(2*dx))*(right - left);
convection = @(right,left) (-c/dx)*(right-left);
convection2 = @(current,left,left_left) (-c/(2*dx)).*(3*current - 4*left + left_left);

%Initial Condition
t = 0;
u = zeros(N,1);
u(1) = sin(A*t); u(N) = 0;
IC = @(t) sin(A*t);



if option == 1
    %% Problem 1
    figure; 
    for i = 1 : nof_t
        
        %-----------------------First Substep---------------------------------%
        t = t + (1/3)*dt;
        f1 = RHS(u(3:N),u(1:N-2));
        u_star(2:N-1,1) = u(2:N-1,1) + a21*dt*f1;
        u_star(1) = IC(t); u_star(N) = 0;
        
        %----------------------Second Substep---------------------------------%
        t = t + (1/3)*dt;
        f2 = RHS(u_star(3:N,1),u_star(1:N-2,1));
        u_star_star(2:N-1,1) = u_star(2:N-1,1) + a32*dt*f2 + zeta2*dt*f1;
        u_star_star(1) = IC(t); u_star_star(N) = 0;
        
        %----------------------Third Substep----------------------------------%
        t = dt*i;
        f3 = RHS(u_star_star(3:N,1), u_star_star(1:N-2,1));
        u(2:N-1) = u_star_star(2:N-1,1) + b3*dt*f3 + zeta3*dt*f2;
        u(1) = IC(t); u(N) = 0;
        
        %     if mod(i,2) == 1
        plot(x,u); title('Pinned Wave'); xlim([0,10]); ylim([0,1.2]);
        pause(0.01)
        %         drawnow
        %     end
    end
    %% Problem 2
    t = 0;
    u = zeros(N,1);
    u(1) = sin(A*t); u(N) = 0;
    figure; 
    for i = 1 : nof_t
        
        %-----------------------First Substep---------------------------------%
        t = t + (1/3)*dt;
        f1 = RHS(u(3:N),u(1:N-2));
        u_star(2:N-1,1) = u(2:N-1,1) + a21*dt*f1;
        u_star(1) = IC(t); u_star(N) = 2*u_star(N-1) - u_star(N-2);
        
        %----------------------Second Substep---------------------------------%
        t = t + (1/3)*dt;
        f2 = RHS(u_star(3:N,1),u_star(1:N-2,1));
        u_star_star(2:N-1,1) = u_star(2:N-1,1) + a32*dt*f2 + zeta2*dt*f1;
        u_star_star(1) = IC(t); u_star_star(N) = 2*u_star_star(N-1) - u_star_star(N-2);
        
        %----------------------Third Substep----------------------------------%
        t = dt*i;
        f3 = RHS(u_star_star(3:N,1), u_star_star(1:N-2,1));
        u(2:N-1) = u_star_star(2:N-1,1) + b3*dt*f3 + zeta3*dt*f2;
        u(1) = IC(t); u(N) = 2*u(N-1) - u(N-2);
        
        %     if mod(i,10) == 1
        plot(x,u); title('Linear Extrapolating'); xlim([0,10]); ylim([0,1.2]);
        pause(0.01)
           
        %         drawnow
        %     end
    end
    %% Problem 3
    t = 0;
    u = zeros(N,1);
    u(1) = sin(A*t); u(N) = 0;
    figure; 
    for i = 1 : nof_t
        
        %-----------------------First Substep---------------------------------%
        t = t + (1/3)*dt;
        f1 = RHS(u(3:N),u(1:N-2));
        u_star(2:N-1,1) = u(2:N-1,1) + a21*dt*f1;
        u_star(1) = IC(t); u_star(N) = 3*u_star(N-1) - 3*u_star(N-2) + u_star(N-3);
        
        %----------------------Second Substep---------------------------------%
        t = t + (1/3)*dt;
        f2 = RHS(u_star(3:N,1),u_star(1:N-2,1));
        u_star_star(2:N-1,1) = u_star(2:N-1,1) + a32*dt*f2 + zeta2*dt*f1;
        u_star_star(1) = IC(t);
        u_star_star(N) = 3*u_star_star(N-1) - 3*u_star_star(N-2) + u_star_star(N-3);
        
        %----------------------Third Substep----------------------------------%
        t = dt*i;
        f3 = RHS(u_star_star(3:N,1), u_star_star(1:N-2,1));
        u(2:N-1) = u_star_star(2:N-1,1) + b3*dt*f3 + zeta3*dt*f2;
        u(1) = IC(t); u(N) = 3*u(N-1) - 3*u(N-2) + u(N-3);
        
        %     if mod(i,2) == 1
        plot(x,u); title('Quadratic Extrapolating'); xlim([0,10]); ylim([0,1.2]);
        pause(0.01)
        %         drawnow
        %     end
    end
    %% Problem 4
    t = 0;
    u = zeros(N,1);
    u(1) = sin(A*t); u(N) = 0;
    figure;
    for i = 1 : nof_t
        
        %-----------------------First Substep---------------------------------%
        t = t + (1/3)*dt;
        f1 = RHS(u(3:N),u(1:N-2));
        u_star(2:N-1,1) = u(2:N-1,1) + a21*dt*f1;
        u_star(1) = IC(t); u_star(N) = u_star(N-1);
        
        %----------------------Second Substep---------------------------------%
        t = t + (1/3)*dt;
        f2 = RHS(u_star(3:N,1),u_star(1:N-2,1));
        u_star_star(2:N-1,1) = u_star(2:N-1,1) + a32*dt*f2 + zeta2*dt*f1;
        u_star_star(1) = IC(t); u_star_star(N) = u_star_star(N-1);
        
        %----------------------Third Substep----------------------------------%
        t = dt*i;
        f3 = RHS(u_star_star(3:N,1), u_star_star(1:N-2,1));
        u(2:N-1) = u_star_star(2:N-1,1) + b3*dt*f3 + zeta3*dt*f2;
        u(1) = IC(t); u(N) = u(N-1);
        
        %     if mod(i,2) == 1
        plot(x,u); ylim([0,1.2]); xlim([0,10]); title('Zero Gradient');
        pause(0.01)
        %         drawnow
        %     end
    end
    %% Problem 5
    t = 0;
    u = zeros(N,1);
    u(1) = sin(A*t); u(N) = 0;
    figure; 
    for i = 1 : nof_t
        
        %-----------------------First Substep---------------------------------%
        t = t + (1/3)*dt;
        f1 = RHS(u(3:N),u(1:N-2));
        u_star(2:N-1,1) = u(2:N-1,1) + a21*dt*f1;
        u_star(1) = IC(t); u_star(N) = -u_star(N-1);
        
        %----------------------Second Substep---------------------------------%
        t = t + (1/3)*dt;
        f2 = RHS(u_star(3:N,1),u_star(1:N-2,1));
        u_star_star(2:N-1,1) = u_star(2:N-1,1) + a32*dt*f2 + zeta2*dt*f1;
        u_star_star(1) = IC(t); u_star_star(N) = -u_star_star(N-1);
        
        %----------------------Third Substep----------------------------------%
        t = dt*i;
        f3 = RHS(u_star_star(3:N,1), u_star_star(1:N-2,1));
        u(2:N-1) = u_star_star(2:N-1,1) + b3*dt*f3 + zeta3*dt*f2;
        u(1) = IC(t); u(N) = -u(N-1);
        
        %     if mod(i,2) == 1
        plot(x,u); title('Antisymmetric'); xlim([0,10]); ylim([0,1.2]);
        pause(0.01)
        %         drawnow
        %     end
    end
    %% Problem 6
    t = 0;
    u = zeros(N,1);
    u(1) = sin(A*t); u(N) = 0;
    figure; 
    for i = 1 : nof_t
        
        %-----------------------First Substep---------------------------------%
        t = t + (1/3)*dt;
        f1 = RHS(u(3:N),u(1:N-2));
        f1_bc = convection(u(N),u(N-1));
        u_star(2:N-1,1) = u(2:N-1,1) + a21*dt*f1;
        u_star(N) = u(N) + a21*dt*f1_bc;
        u_star(1) = IC(t);
        
        %----------------------Second Substep---------------------------------%
        t = t + (1/3)*dt;
        f2 = RHS(u_star(3:N,1),u_star(1:N-2,1));
        f2_bc = convection(u_star(N),u_star(N-1));
        u_star_star(2:N-1,1) = u_star(2:N-1,1) + a32*dt*f2 + zeta2*dt*f1;
        u_star_star(N) = u_star(N) + a32*dt*f2_bc + zeta2*dt*f1_bc;
        u_star_star(1) = IC(t);
        
        
        %----------------------Third Substep----------------------------------%
        t = dt*i;
        f3 = RHS(u_star_star(3:N,1), u_star_star(1:N-2,1));
        f3_bc = convection(u_star_star(N),u_star_star(N-1));
        u(2:N-1) = u_star_star(2:N-1,1) + b3*dt*f3 + zeta3*dt*f2;
        u(N) = u_star_star(N) + b3*dt*f3_bc + zeta3*dt*f2_bc;
        u(1) = IC(t);
        
        %     if mod(i,2) == 1
        plot(x,u); xlim([0,10]); ylim([0,1.2]);
        title('First Order Upwind');
        pause(0.01)
        %         drawnow
        %     end
    end
    %% Problem 7
    t = 0;
    u = zeros(N,1);
    u(1) = sin(A*t); u(N) = 0;
    figure; 
    for i = 1 : nof_t
        
        %-----------------------First Substep---------------------------------%
        t = t + (1/3)*dt;
        f1 = RHS(u(3:N),u(1:N-2));
        f1_bc = convection2(u(N),u(N-1),u(N-2));
        u_star(2:N-1,1) = u(2:N-1,1) + a21*dt*f1;
        u_star(N) = u(N) + a21*dt*f1_bc;
        u_star(1) = IC(t);
        
        %----------------------Second Substep---------------------------------%
        t = t + (1/3)*dt;
        f2 = RHS(u_star(3:N,1),u_star(1:N-2,1));
        f2_bc = convection2(u_star(N),u_star(N-1),u_star(N-1));
        u_star_star(2:N-1,1) = u_star(2:N-1,1) + a32*dt*f2 + zeta2*dt*f1;
        u_star_star(N) = u_star(N) + a32*dt*f2_bc + zeta2*dt*f1_bc;
        u_star_star(1) = IC(t);
        
        
        %----------------------Third Substep----------------------------------%
        t = dt*i;
        f3 = RHS(u_star_star(3:N,1), u_star_star(1:N-2,1));
        f3_bc = convection2(u_star_star(N),u_star_star(N-1),u_star_star(N-2));
        u(2:N-1) = u_star_star(2:N-1,1) + b3*dt*f3 + zeta3*dt*f2;
        u(N) = u_star_star(N) + b3*dt*f3_bc + zeta3*dt*f2_bc;
        u(1) = IC(t);
        
        %     if mod(i,10) == 1
        plot(x,u); title('Second Order Upwind'); xlim([0,10]); ylim([0,1.2]);
        pause(0.01)
        %         drawnow
        %     end
    end
else
    if problem == 1
        %%
        figure
        for i = 1 : nof_t
            
            %-----------------------First Substep---------------------------------%
            t = t + (1/3)*dt;
            f1 = RHS(u(3:N),u(1:N-2));
            u_star(2:N-1,1) = u(2:N-1,1) + a21*dt*f1;
            u_star(1) = IC(t); u_star(N) = 0;
            
            %----------------------Second Substep---------------------------------%
            t = t + (1/3)*dt;
            f2 = RHS(u_star(3:N,1),u_star(1:N-2,1));
            u_star_star(2:N-1,1) = u_star(2:N-1,1) + a32*dt*f2 + zeta2*dt*f1;
            u_star_star(1) = IC(t); u_star_star(N) = 0;
            
            %----------------------Third Substep----------------------------------%
            t = dt*i;
            f3 = RHS(u_star_star(3:N,1), u_star_star(1:N-2,1));
            u(2:N-1) = u_star_star(2:N-1,1) + b3*dt*f3 + zeta3*dt*f2;
            u(1) = IC(t); u(N) = 0;
            
            %     if mod(i,2) == 1
            plot(x,u); ylim([0,1.2]); xlim([0,10]); title('Pinned Wave');
            pause(0.01)
            %         drawnow
            %     end
        end
        figure
    elseif problem == 2
        %%
        figure;
        for i = 1 : nof_t
            
            %-----------------------First Substep---------------------------------%
            t = t + (1/3)*dt;
            f1 = RHS(u(3:N),u(1:N-2));
            u_star(2:N-1,1) = u(2:N-1,1) + a21*dt*f1;
            u_star(1) = IC(t); u_star(N) = 2*u_star(N-1) - u_star(N-2);
            
            %----------------------Second Substep---------------------------------%
            t = t + (1/3)*dt;
            f2 = RHS(u_star(3:N,1),u_star(1:N-2,1));
            u_star_star(2:N-1,1) = u_star(2:N-1,1) + a32*dt*f2 + zeta2*dt*f1;
            u_star_star(1) = IC(t); u_star_star(N) = 2*u_star_star(N-1) - u_star_star(N-2);
            
            %----------------------Third Substep----------------------------------%
            t = dt*i;
            f3 = RHS(u_star_star(3:N,1), u_star_star(1:N-2,1));
            u(2:N-1) = u_star_star(2:N-1,1) + b3*dt*f3 + zeta3*dt*f2;
            u(1) = IC(t); u(N) = 2*u(N-1) - u(N-2);
            
            %     if mod(i,10) == 1
            plot(x,u); ylim([0,1.2]); xlim([0,10]); title('Linear Extrapolating');
            pause(0.01)
            %         drawnow
            %     end
        end

    elseif problem ==3
        %%
        figure
        for i = 1 : nof_t
            
            %-----------------------First Substep---------------------------------%
            t = t + (1/3)*dt;
            f1 = RHS(u(3:N),u(1:N-2));
            u_star(2:N-1,1) = u(2:N-1,1) + a21*dt*f1;
            u_star(1) = IC(t); u_star(N) = 3*u_star(N-1) - 3*u_star(N-2) + u_star(N-3);
            
            %----------------------Second Substep---------------------------------%
            t = t + (1/3)*dt;
            f2 = RHS(u_star(3:N,1),u_star(1:N-2,1));
            u_star_star(2:N-1,1) = u_star(2:N-1,1) + a32*dt*f2 + zeta2*dt*f1;
            u_star_star(1) = IC(t);
            u_star_star(N) = 3*u_star_star(N-1) - 3*u_star_star(N-2) + u_star_star(N-3);
            
            %----------------------Third Substep----------------------------------%
            t = dt*i;
            f3 = RHS(u_star_star(3:N,1), u_star_star(1:N-2,1));
            u(2:N-1) = u_star_star(2:N-1,1) + b3*dt*f3 + zeta3*dt*f2;
            u(1) = IC(t); u(N) = 3*u(N-1) - 3*u(N-2) + u(N-3);
            
            %     if mod(i,2) == 1
            plot(x,u); ylim([0,1.2]); xlim([0,10]); title('Quadratic Extrapolating');
            pause(0.01)
            %         drawnow
            %     end
        end
        
    elseif problem == 4
        %%
        figure; 
        for i = 1 : nof_t
            
            %-----------------------First Substep---------------------------------%
            t = t + (1/3)*dt;
            f1 = RHS(u(3:N),u(1:N-2));
            u_star(2:N-1,1) = u(2:N-1,1) + a21*dt*f1;
            u_star(1) = IC(t); u_star(N) = u_star(N-1);
            
            %----------------------Second Substep---------------------------------%
            t = t + (1/3)*dt;
            f2 = RHS(u_star(3:N,1),u_star(1:N-2,1));
            u_star_star(2:N-1,1) = u_star(2:N-1,1) + a32*dt*f2 + zeta2*dt*f1;
            u_star_star(1) = IC(t); u_star_star(N) = u_star_star(N-1);
            
            %----------------------Third Substep----------------------------------%
            t = dt*i;
            f3 = RHS(u_star_star(3:N,1), u_star_star(1:N-2,1));
            u(2:N-1) = u_star_star(2:N-1,1) + b3*dt*f3 + zeta3*dt*f2;
            u(1) = IC(t); u(N) = u(N-1);
            
            %     if mod(i,2) == 1
            plot(x,u); ylim([0,1.2]); xlim([0,10]); title('Zero Gradient');
            pause(0.01)
            %         drawnow
            %     end
        end
        
    elseif problem == 5
        %%
        figure; %Problem5
        for i = 1 : nof_t
            
            %-----------------------First Substep---------------------------------%
            t = t + (1/3)*dt;
            f1 = RHS(u(3:N),u(1:N-2));
            u_star(2:N-1,1) = u(2:N-1,1) + a21*dt*f1;
            u_star(1) = IC(t); u_star(N) = -u_star(N-1);
            
            %----------------------Second Substep---------------------------------%
            t = t + (1/3)*dt;
            f2 = RHS(u_star(3:N,1),u_star(1:N-2,1));
            u_star_star(2:N-1,1) = u_star(2:N-1,1) + a32*dt*f2 + zeta2*dt*f1;
            u_star_star(1) = IC(t); u_star_star(N) = -u_star_star(N-1);
            
            %----------------------Third Substep----------------------------------%
            t = dt*i;
            f3 = RHS(u_star_star(3:N,1), u_star_star(1:N-2,1));
            u(2:N-1) = u_star_star(2:N-1,1) + b3*dt*f3 + zeta3*dt*f2;
            u(1) = IC(t); u(N) = -u(N-1);
            
            %     if mod(i,2) == 1
            plot(x,u); ylim([0,1.2]); xlim([0,10]); title('Antisymmetric');
            pause(0.01)
            %         drawnow
            %     end
        end
    elseif problem == 6
        %%
        figure; %Problem6
        for i = 1 : nof_t
            
            %-----------------------First Substep---------------------------------%
            t = t + (1/3)*dt;
            f1 = RHS(u(3:N),u(1:N-2));
            f1_bc = convection(u(N),u(N-1));
            u_star(2:N-1,1) = u(2:N-1,1) + a21*dt*f1;
            u_star(N) = u(N) + a21*dt*f1_bc;
            u_star(1) = IC(t);
            
            %----------------------Second Substep---------------------------------%
            t = t + (1/3)*dt;
            f2 = RHS(u_star(3:N,1),u_star(1:N-2,1));
            f2_bc = convection(u_star(N),u_star(N-1));
            u_star_star(2:N-1,1) = u_star(2:N-1,1) + a32*dt*f2 + zeta2*dt*f1;
            u_star_star(N) = u_star(N) + a32*dt*f2_bc + zeta2*dt*f1_bc;
            u_star_star(1) = IC(t);
            
            
            %----------------------Third Substep----------------------------------%
            t = dt*i;
            f3 = RHS(u_star_star(3:N,1), u_star_star(1:N-2,1));
            f3_bc = convection(u_star_star(N),u_star_star(N-1));
            u(2:N-1) = u_star_star(2:N-1,1) + b3*dt*f3 + zeta3*dt*f2;
            u(N) = u_star_star(N) + b3*dt*f3_bc + zeta3*dt*f2_bc;
            u(1) = IC(t);
            
            %     if mod(i,2) == 1
            plot(x,u);
            ylim([0,1.2]); xlim([0,10]); title('First Order Upwind');
            pause(0.01)
            %         drawnow
            %     end
        end
    elseif problem == 7 
        %%
        figure; %Problem7
        for i = 1 : nof_t
            
            %-----------------------First Substep---------------------------------%
            t = t + (1/3)*dt;
            f1 = RHS(u(3:N),u(1:N-2));
            f1_bc = convection2(u(N),u(N-1),u(N-2));
            u_star(2:N-1,1) = u(2:N-1,1) + a21*dt*f1;
            u_star(N) = u(N) + a21*dt*f1_bc;
            u_star(1) = IC(t);
            
            %----------------------Second Substep---------------------------------%
            t = t + (1/3)*dt;
            f2 = RHS(u_star(3:N,1),u_star(1:N-2,1));
            f2_bc = convection2(u_star(N),u_star(N-1),u_star(N-1));
            u_star_star(2:N-1,1) = u_star(2:N-1,1) + a32*dt*f2 + zeta2*dt*f1;
            u_star_star(N) = u_star(N) + a32*dt*f2_bc + zeta2*dt*f1_bc;
            u_star_star(1) = IC(t);
            
            
            %----------------------Third Substep----------------------------------%
            t = dt*i;
            f3 = RHS(u_star_star(3:N,1), u_star_star(1:N-2,1));
            f3_bc = convection2(u_star_star(N),u_star_star(N-1),u_star_star(N-2));
            u(2:N-1) = u_star_star(2:N-1,1) + b3*dt*f3 + zeta3*dt*f2;
            u(N) = u_star_star(N) + b3*dt*f3_bc + zeta3*dt*f2_bc;
            u(1) = IC(t);
            
            %     if mod(i,10) == 1
            plot(x,u); ylim([0,1.2]); xlim([0,10]); title('Second Order Upwind');
            pause(0.01)
            %         drawnow
            %     end
        end
    end
end








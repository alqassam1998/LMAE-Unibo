%% SYSTEM THEORY AND ADVANCED CONTROL - Academic Year 2024/2025
% EXERCISE SESSION V
% Reachability Discrete Time

clear all; clc; close all;

%% System Model

% Details about the state variables and control inputs meaning in the
% addendum PDF

A = [0.8, -0.2, -0.1;
    0.1, 0.5, 0.1;
    0.1, - 0.1, 0.5
    ];

B = [2, 0, -0.1, 0;
    0, 0.3, -0.4, - 0.5;
    0.5, 0.2, -0.5, -0.4
    ];

eigenvalues = eigs(A);
disp('Eigenvalues of A')
disp(eigenvalues)
disp('Rank of B')
disp(rank(B))
%% Computation of the Reachability/Controllability Matrix
R = [B A*B A*A*B]; 
R_matlab = ctrb(A, B);

disp('Reachability Matrix by Hand')
disp(R)
disp('Rank of R:')
disp(rank(R))

disp('Reachability Matrix via MATLAB command')
disp(R_matlab)
disp('Rank of R matlab command:')
disp(rank(R_matlab))

%% System behavior with no input u = 0

x0 = [200;
    5; 
    40
    ];

n_steps = 20;
xx = zeros(size(x0, 1), n_steps);
xx(:, 1) = x0;

for k=2:(n_steps + 1)
    xx(:, k) = A * xx(:, k - 1);
end

figure;
ax1 = gca;
hold(ax1, 'all');
stairs(0:n_steps', xx', LineWidth=1.5);
legend('Chickens, $x_1$', 'Wolves, $x_2$', 'Boars, $x_3$', interpreter='latex', ...
    location='northeast')
grid();
xlabel('Steps')
ylabel('System States')

%% Computation of the Control Input u

x0 = [200;
    5; 
    40
    ];

x_ref = [350;
    15;
    50
    ];

k_ref = 3;

n = size(A, 1);
m = size(B, 2);

R_ref = zeros(n, n*m);

for k=1:k_ref
    R_ref(:, m*(k-1) + 1:m*k) = A^(k - 1)*B;
end

uu = R_ref' * pinv(R_ref * R_ref')*(x_ref - A^k_ref * x0);
% disp(uu)
uu_seq = fliplr(reshape(uu, m, k_ref));
% disp(uu_seq)


%% Simulate the Results
close all;
n_steps = k_ref;

xx = zeros(size(x0, 1) ,n_steps + 1);
xx(:,1) = x0;

for k=2:n_steps + 1
   xx(:,k) = A*xx(:, k-1)+ B*uu_seq(:, k-1);
end

figure;
ax2 = gca;
hold(ax2, 'all');
stairs(ax2, 0:n_steps', xx', LineWidth=2);
plot(ax2, 0:n_steps, x_ref .* ones(size(x_ref, 1), n_steps + 1), LineStyle='--', LineWidth=1.5);
legend(ax2, 'Chickens, $x_1$', 'Wolves, $x_2$', 'Boars, $x_3$',...
    '$x_{ref,1}$', '$x_{ref,2}$', '$x_{ref,3}$', interpreter='latex', ...
    location='northwest')
grid(ax2);
xlabel(ax2, 'Steps')
ylabel(ax2, 'System States')
title(ax2, 'Reach of the target with arbitrary k')



%% Can we lower the k_ref?
% What is the system Reachability Index?

rank(B)

r_index = 1;

k_ref = r_index;

% Let's compute the new control input 

x0 = [200;
    5; 
    40
    ];

x_ref = [500;
    150;
    200
    ];

R_ref = B;

uu = R_ref' * pinv(R_ref * R_ref')*(x_ref - A^k_ref * x0);
% disp(uu)
uu_seq = fliplr(reshape(uu, m, k_ref));
% disp(uu_seq)

%% Simulate the Results
n_steps = r_index;

xx = zeros(size(x0, 1) ,n_steps + 1);
xx(:,1) = x0;

for k=2:n_steps + 1
   xx(:,k) = A*xx(:, k-1)+ B*uu_seq(:, k-1);
end

figure;
ax3 = gca;
hold(ax3, 'all');
stairs(ax3, 0:n_steps', xx', LineWidth=1.5);
plot(ax3, 0:n_steps, x_ref .* ones(size(x_ref, 1), n_steps + 1), LineStyle='--', LineWidth=1.5);
legend(ax3, 'Chickens, $x_1$', 'Wolves, $x_2$', 'Boars, $x_3$',...
    '$x_{ref,1}$', '$x_{ref,2}$', '$x_{ref,3}$', interpreter='latex', ...
    location='northwest')
grid(ax3);
xlabel(ax3, 'Steps')
ylabel(ax3, 'System States')
title(ax3, 'Reach of the target state with Reach. Index')

%% Is it possible with a fewer number of inputs
% to reach the exact reference state in one step?
% No, it is not. We know that if the rank of B < n, the best we can get is
% to reach the state after n steps (if the system is C.R.)

% Suppose to reduce the number of inputs from 4 to 2
clc; clear all; 

A = [0.8, -0.2, -0.1;
    0.1, 0.5, 0.1;
    0.1, - 0.1, 0.5
    ];

B = [2, 0;
    0, 0.3;
    0.5, 0.2
    ];

disp('Rank of the matrix B');
rank(B)
% Let's check again the controllability matrix
R = ctrb(A, B);
disp('Controllability Matrix rank');
rank(R)

R_manual = [B A*B];
disp('Controllability Matrix rank (by hand)')
rank(R_manual)

r_index = 1;
k_ref = 3; % Is it possible if rank B < n?
m = size(B, 2);

% Let's compute the new control input 
x0 = [200;
    5; 
    40
    ];

x_ref = [350;
    15;
    50
    ];

n = size(A, 1);
m = size(B, 2);

R_ref = zeros(n, n*m);

% Populate the Reachability matrix
for k=1:k_ref
    R_ref(:, m*(k-1) + 1:m*k) = A^(k - 1)*B;
end

% Computation of the Control Input u
uu = R_ref' * pinv(R_ref * R_ref')*(x_ref - A^k_ref * x0);
disp(uu)

% Rearrange the vector to obtain the correct dimension
uu_seq = fliplr(reshape(uu, m, k_ref));
% disp(uu_seq)

n_steps = k_ref;

xx = zeros(size(x0, 1) ,n_steps + 1);
xx(:,1) = x0;

for k=2:n_steps + 1
   xx(:,k) = A*xx(:, k-1)+ B*uu_seq(:, k-1);
end

figure;
ax4 = gca;
hold(ax4, 'all');
stairs(ax4, 0:n_steps', xx', LineWidth=1.5);
plot(ax4, 0:n_steps, x_ref .* ones(size(x_ref, 1), n_steps + 1), LineStyle='--', LineWidth=1.5);
legend(ax4, 'Chickens, $x_1$', 'Wolves, $x_2$', 'Boars, $x_3$',...
    '$x_{ref,1}$', '$x_{ref,2}$', '$x_{ref,3}$', interpreter='latex', ...
    location='northwest')
grid(ax4);
xlabel(ax4, 'Steps')
ylabel(ax4, 'System States')
title(ax4, 'Reach of the target state with fewer inputs')

%% Extra - Add an Arbitrary vector into the Kernel of R_ref
% You can add any other element in the kernel of R_ref without affecting
% the end state but perhaps affecting the trajectory and the control energy

% 1) Compare uu_seq computed before with another uu_seq + deltaU, where deltaU
% is an arbitrary vector from the kernel of R_ref

% hint 
%deltaU = ker_Rref * scalar * rand
% u* + deltaU

% 2) Plot and compare the two norms

%% Extra - Compute the Equilibrium input U_eq for X_ref
% The computed control input u, is to 'reach' the target state, can we
% actually hold the state at x_reference as t -> infinity?

% u* reachability control input
% u_eq keep the state to x_ref

% 1) Check if the computed uu_seq is able to maintain the system at x_ref

% hint - create another vector with uu_seq and simulate for k > k_ref

% 2) Compute the equilibrium control input to reach x_ref
% hint - rely on the eq. pair (x_eq, u_eq) from the theory of DT systems.

%% Extra - Exploring Redundancy.

% What about, for example, lowering the amount of poison realeased in the
% environment?

% 1) With the same idea exploited before (thanks to the kernel of R_ref) try
% to compute a control input u reaching the same x_ref but lowering the
% amount of poison released (i.e. input u3).

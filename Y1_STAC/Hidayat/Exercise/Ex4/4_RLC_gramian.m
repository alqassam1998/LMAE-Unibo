%% SYSTEM THEORY AND ADVANCED CONTROL - Academic Year 2024/2025
% EXERCISE SESSION V
% Reachability Continuous Time
clc; clear all; close all;

L=1;
C=1e-1;
R = sqrt(L/C)*1;   %try this later
% R = 10;

% x_1 = i_L
% x_2 = v_c

A = [-(2*R)/L, 1/L
    - 1/C, 0];

B = [R/L;
    1/C
    ];

C = [1, 0];


[V, lambda] = eig(A);

% A_tilde = pinv(V) * A * V

kerA = null(lambda(1) * eye(2) - A);
disp('Kernel of A')
disp(kerA);

disp('Eigenvalues of A');
disp(lambda);

disp('Eigenvectors of A')
disp(V)


% sys = ss(A, B, eye(size(A, 1)), zeros(size(A, 1), size(B, 2)));
% bode(sys);
% Gs = tf(sys);
% Gs


disp(lambda);
R=ctrb(A,B);
disp('Rank of R')
disp(rank(R))

%% System behavior in absence of inputs

x0 = [0.5; 
    0.5];

t_end = 10;

model = @(t, x) A * x;

dt = 1e-4;
[time, xx] = ode45(model, 0:dt:t_end, x0);

figure;
ax0 = gca;
hold(ax0, "on")
plot(ax0, time, xx(:, 1), LineWidth=1.5, DisplayName='$x_1$')
plot(ax0, time, xx(:, 2), LineWidth=1.5, DisplayName='$x_2$')
ylabel(ax0, 'System States', Interpreter='latex')
xlabel(ax0, 'Time (s)', Interpreter='latex')
grid(ax0)
legend(ax0, Interpreter='latex')
%%
close all;
% control goal:
Tf = 1;  %decrease this!

x_ref = [10;
    15];

%computation of the Gramian at Tf
F = @(t,x) reshape(expm(A*(Tf-t))*B*B'*expm(A'*(Tf-t)), 4, 1);

time_span = 0:Tf*dt:Tf;

[~, w] = ode45(F, time_span, zeros(4, 1));

W = reshape(w(end,:)', 2, 2);
disp('Determinant of W')
det(W)
invW = pinv(W);

%construction of u
u = @(t, x0) B'* expm(A' * (Tf-t)) * invW * (x_ref - expm(A*Tf) * x0);

F = @(t,x) A*x+B*u(t,x0);

% model = @(t, x, x0) [ % u ....
%     %F...];
[t,x] = ode45(F, time_span, x0);

U = zeros(size(t));
for i=1:length(t)
    U(i)=u(t(i),x0);
end

% First subplot
figure;
ax1 = subplot(2,1,1); 
hold(ax1, 'on');
plot(ax1, t, x(:, 1), 'linewidth', 1.5, DisplayName='$x_1$');
plot(ax1, t, x(:, 2), 'linewidth', 1.5, DisplayName='$x_2$');
plot(ax1, t, x_ref(1) * ones(size(t)), LineStyle="--", DisplayName='$x_{1, ref}$')
plot(ax1, t, x_ref(2) * ones(size(t)), LineStyle="--", DisplayName='$x_{2, ref}$')
xlabel(ax1, 'Time (s)', interpreter='latex')
ylabel(ax1, '$x(t)$', Interpreter='latex');
legend(ax1, interpreter='latex', Location='southwest')
grid(ax1)

% Second subplot
ax2 = subplot(2,1,2);  
plot(ax2, t, U, 'linewidth', 1.5, DisplayName='$u(t)$');
hold(ax2, 'on');
hold(ax2, 'off');
xlabel(ax2, 'Time (s)', Interpreter='latex');
ylabel(ax2, '$u(t)$', Interpreter='latex');
legend(ax2, Interpreter="latex")
grid(ax2)

figure;
ax4 = gca;
hold(ax4, "on");
plot(ax4, x(:, 1), x(:, 2), LineWidth=1.5);
xlabel(ax4, '$x_1$', interpreter='latex');
ylabel(ax4, '$x_2$', Interpreter='latex');
grid(ax4)
title(ax4, 'Phase Space Plot', Interpreter='latex')


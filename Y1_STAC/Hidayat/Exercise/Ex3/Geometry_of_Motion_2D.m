%% SYSTEM THEORY AND ADVANCED CONTROL - Academic Year 2024/2025
% EXERCISE SESSION III 
% Geometry of Motion - 2D Evolution
clc; clear variables; close all;

% 1) System with Two Stable Eigenvalues (Stable Nodes)
% Define system dynamics with Hurwitz matrix
disp('1) Two stable eigenvalues: Stable nodes');

A = [-2.2 -1.0198; -0.9414 -1.8]; % Example matrix A with stable eigenvalues
[V, D] = eig(A); % Calculate eigenvalues and eigenvectors

disp('Eigenvalues of A:');
disp(D);

% Define grid range for state-space
L = 8;
[x1, x2] = meshgrid(linspace(-L, L, 40)); % Create grid for phase portrait
x1_dot = A(1,1)*x1 + A(1,2)*x2;
x2_dot = A(2,1)*x1 + A(2,2)*x2;

% Plot phase portrait
figure;
ax1 = gca;
hold(ax1, 'all');
quiver(ax1, x1, x2, x1_dot, x2_dot, 'AutoScaleFactor', 0.7, 'LineWidth', 1.5, 'Color', [0.3, 0.6, 0.8], DisplayName='Vector Field');
axis([-L L -L L]);
title(ax1, 'Phase Portrait with Stable Eigenvalues');
xlabel(ax1, 'x1'); 
ylabel(ax1, 'x2');

% Plot eigenspaces
l = linspace(-L, L, 10);
v1 = (V(:,1)./norm(V(:,1))) * l;
v2 = (V(:,2)./norm(V(:,2))) * l;
plot(ax1, v1(1,:), v1(2,:), 'LineWidth', 2, DisplayName='Stable Node Eigenvector 1');
plot(ax1, v2(1,:), v2(2,:), 'LineWidth', 2, DisplayName='Stable Node Eigenvector 2');
legend(ax1, 'Interpreter', 'latex', Location='southeastoutside')
grid on;

%% 1.1) Animate Solution with Initial Condition
x0 = -5*V(:,1)./norm(V(:,1)) + -5*V(:,2)./norm(V(:,2)); % Initial condition

eps = 0.01;
Tf = log(norm(x0)/eps) / -max(eig(A));
t = 0:0.08:Tf;

p = [];
q = [];
% Plot animation of solution
initial_point = plot(ax1, x0(1), x0(2), 'ko', 'MarkerFaceColor', 'k', DisplayName='$x_0$');
for i = 2:length(t)
    x = solveLinearSystem(A, t(1:i), x0); % Solve system from t=0 to t(i)
    if i > 2 
        delete(p); 
        delete(q); 
    end
    p = plot(ax1, x(1,:), x(2,:),'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'r', DisplayName='Trajectory');
    q = plot(ax1, x(1,end), x(2,end), 'bo', 'MarkerFaceColor', 'b', DisplayName='$x_{\infty}$');
    pause(0.05);
end

%% 1.2) Plot Family of Solutions

% Clear previous plots
delete(p); 
delete(q); 
delete(initial_point);

theta = linspace(0, 2*pi, 20); % Parametric family of initial conditions
for i = 1:length(theta)
    x0 = 7 * [cos(theta(i)); sin(theta(i))];
    x = solveLinearSystem(A, t, x0);
    plot(ax1, x(1,1), x0(2,1), 'bo', 'MarkerFaceColor', 'b', DisplayName='Initial Point'); % Initial point
    plot(ax1, x(1,:), x(2,:), 'LineWidth', 1.5, LineStyle='--');
    delete(legend)
end
legend(ax1, '', 'Stable Node Eigenvector 1', 'Stable Node Eigenvector 2');
title(ax1, 'Family of Solutions for Stable Node');

%% 2) System with One Stable and One Unstable Eigenvalue (Saddle Point)
clc; clear variables; close all;

A = [-0.2 -1.8; -1.2 -0.8];
[V, D] = eig(A); % Compute eigenvalues and eigenvectors
disp('2) Saddle Point: Eigenvalues of A');
disp(D);

L = 8; % Range limit
[x1, x2] = meshgrid(linspace(-L, L, 40));
x1_dot = A(1,1)*x1 + A(1,2)*x2;
x2_dot = A(2,1)*x1 + A(2,2)*x2;

% Phase portrait for saddle point
figure;
ax2 = gca;
hold(ax2, 'all');
quiv = quiver(ax2, x1, x2, x1_dot, x2_dot, 'AutoScaleFactor', 1, 'LineWidth', 1.5, 'Color', [0.3, 0.5, 0.7], DisplayName='Vector Field');
axis([-L L -L L]);
title(ax2, 'Phase Portrait with Saddle Point');
xlabel(ax2, 'x1'); ylabel(ax2, 'x2');

l = linspace(-L, L, 10); % Range for eigenspaces
v1 = (V(:,1)./norm(V(:,1))) * l;
v2 = (V(:,2)./norm(V(:,2))) * l;

plot(ax2, v1(1,:), v1(2,:), 'LineWidth', 2, DisplayName='Unstable Eigenvector');
plot(ax2, v2(1,:), v2(2,:), 'LineWidth', 2, DisplayName='Stable Eigenvector');
legend(ax2, 'interpreter', 'latex', Location='southeastoutside');
grid on;

%% 2.1) Animate Solution with Initial Condition for Saddle
x0 = -7*V(:,2) + 0.6*V(:,1);
eps = 1/50;
Tf = max(log(norm(x0)/eps)/-min(D(D<0)), log(norm(x0)*sqrt(2)*L)/max(D(D>0)));
t = 0:0.08:Tf;

p = [];
q = [];
initial_point = plot(ax2, x0(1), x0(2), 'ko', 'MarkerFaceColor', 'k', DisplayName='$x_0$');
for i = 2:length(t)
    x = solveLinearSystem(A, t(1:i), x0);
    if i > 2 
        delete(p); 
        delete(q); 
    end
    p = plot(ax2, x(1,:), x(2,:), 'LineWidth', 1.5,LineStyle='--', DisplayName='Trajectory', Color='red');
    q = plot(ax2, x(1,end), x(2,end), 'bo', 'MarkerFaceColor', 'b', DisplayName='$x_{\infty}$');
    pause(0.05);
    
end

%% 2.2) Plot Family of Solutions for Saddle Point

% Delete previous Elements on the plot
delete(p); 
delete(q);
delete(initial_point);
delete(legend);
delete(quiv);

theta = linspace(0, 2*pi, 20); % parametrized initial conditions
for i = 1:length(theta)
    x0 = 5 * [cos(theta(i)); sin(theta(i))];
    x = solveLinearSystem(A, t, x0);
    plot(ax2, x(1,1), x0(2,1), 'bo', 'MarkerFaceColor', 'b', DisplayName='$x_0$'); % Initial point
    plot(ax2, x(1,:), x(2,:), 'LineWidth', 2, LineStyle='--', DisplayName='');
end
legend(ax2, 'Unstable Eigenvector', 'Stable Eigenvector')
title(ax2, 'Family of Solutions for Saddle Point');

%% 3) 2 Complex Eigenvalues
clear; clc; close all;

alpha = - 0.5;  % Complex number z = alpha + omega
omega = 3;

% Let's construct a matrix having complex conjugate eigenvalues
V = [-1, 3;
    1, 3/2];

S = [alpha, omega; 
    - omega, alpha];

A = V * S * pinv(V);

[V_comp, D] = eig(A);

disp('Eigenvalues of A:');
disp(D);

L = 8; % Range limit
[x1, x2] = meshgrid(linspace(-L, L, 40));
x1_dot = A(1,1)*x1 + A(1,2)*x2;
x2_dot = A(2,1)*x1 + A(2,2)*x2;

% Phase portrait for complex eigenvalues
figure;
ax3 = gca;
hold(ax3, 'all');
quiv = quiver(ax3, x1, x2, x1_dot, x2_dot, 'AutoScaleFactor', 1, 'LineWidth', 1.5, 'Color', [0.3, 0.5, 0.7], DisplayName='Vector Field');
axis([-L L -L L]);
title(ax3, 'Phase Portrait with Complex Eigenvalues');
xlabel(ax3, 'x1'); 
ylabel(ax3, 'x2');

l = linspace(-L, L, 10); % Range for eigenspaces
v1 = (V(:,1)./norm(V(:,1))) * l;
v2 = (V(:,2)./norm(V(:,2))) * l;

plot(ax3, v1(1,:), v1(2,:), 'LineWidth', 2, DisplayName='$v_1$');
plot(ax3, v2(1,:), v2(2,:), 'LineWidth', 2, DisplayName='$v_2$');
legend(ax3, 'interpreter', 'latex', Location='southeastoutside');
grid on;

%% 3.1) Animate the solution
x0 = - 2*V(:,2)./norm(V(:,2))+ 2*V(:,1)./norm(V(:,1));

eps = 1/50;
Tf = 1.2*max(log(norm(x0)/eps)./alpha, log(norm(x0)*sqrt(2)*L)./-alpha);
t = 0:0.08:Tf;

if Tf==Inf
    Tf=4;
end

p = [];
q = [];
initial_point = plot(ax3, x0(1), x0(2), 'ko', 'MarkerFaceColor', 'k', DisplayName='$x_0$');
for i = 2:length(t)
    x = solveLinearSystem(A, t(1:i), x0);
    if i > 2 
        delete(p); 
        delete(q); 
    end
    p = plot(ax3, x(1,:), x(2,:), 'LineWidth', 1.5,LineStyle='--', DisplayName='Trajectory', Color='red');
    q = plot(ax3, x(1,end), x(2,end), 'bo', 'MarkerFaceColor', 'b', DisplayName='$x_{\infty}$');
    pause(0.02);
end

if alpha<0
    title(ax3, 'Stable focus');
else
    title(ax3, 'Unstable focus');
end

%% 2.2) Plot Family of Solutions for Complex Eigenvalues
% Delete previous Elements on the plot
delete(p); 
delete(q);
delete(initial_point);
delete(legend);
delete(quiv);

theta = linspace(0, 2*pi, 10); % parametrized initial conditions
for i = 1:length(theta)
    x0 = 5 * [cos(theta(i)); sin(theta(i))];
    x = solveLinearSystem(A, t, x0);
    plot(ax3, x(1,1), x0(2,1), 'bo', 'MarkerFaceColor', 'b', DisplayName='$x_0$'); % Initial point
    plot(ax3, x(1,:), x(2,:), 'LineWidth', 2, LineStyle='--', DisplayName='');
end
title(ax3, 'Family of Solutions for Complex Eigenvalues');

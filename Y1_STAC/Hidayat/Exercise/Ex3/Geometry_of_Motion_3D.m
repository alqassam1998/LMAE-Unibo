%% SYSTEM THEORY AND ADVANCED CONTROL - Academic Year 2024/2025
% EXERCISE SESSION III 
% Geometry of Motion - 3D Evolution

% 1) 3D with Invariant Eigenspaces
clear; close all; clc;

A = [-1.7986 -.9793 0;
     -.9797 -2.2014 0;
    0 0 -10];

[V,D] = eig(A);

disp('Eigenvalues of A:');
disp(D);

% plot the eigenvectors (exactly as before)
L = 8;
l = linspace(-L,L,10).*sqrt(2); %range for the eigenspaces
v1 = ( V(:,1)./norm(V(:,1)) ) * l; % normalized eigenvector times l
v2 = ( V(:,2)./norm(V(:,2)) ) * l;
v3 = ( V(:,3)./norm(V(:,3)) ) * l;

figure;
ax3 = gca;
plot3(ax3, v1(1,:), v1(2,:), v1(3,:), ...
        v2(1,:), v2(2,:), v2(3,:), ...
        v3(1,:), v3(2,:), v3(3,:), 'linewidth', 2);
legend(ax3, 'v1', 'v2', 'v3', 'Interpreter', 'Latex');
xlabel('x1');
ylabel('x2');
zlabel('x3');
grid(ax3);
axis(1.1 * [-L L -L L -L L]);


% Color the eigenplane spanned by v1 and v2, without adding it to the legend
patch(ax3, [v1(1,1) v2(1,1) v1(1,end) v2(1,end)], ...
           [v1(2,1) v2(2,1) v1(2,end) v2(2,end)], ...
           [v1(3,1) v2(3,1) v1(3,end) v2(3,end)], 'b', ...
           'FaceAlpha', 0.2, 'HandleVisibility', 'off');

%% 1b) Animate the 3D plot for Stable Invariants
x0 = 4 * V(:,1) + 4 * V(:,2) + 4 * V(:,3);  % initial condition
eps = 1/50;
Tf = abs(log(norm(x0) / eps) / -max(D(D ~= 0))); % where D are eigenvalues
t = 0:0.05:Tf;  % time span

hold on;
delete(legend);
% Plot the initial point
initial_point = plot3(ax3, x0(1), x0(2), x0(3), 'Marker','o', 'Color', 'black', 'DisplayName', '$x_0$');
p = [];
q = [];
for i = 2:length(t)
    x = solveLinearSystem(A, t(1:i), x0);
    if i > 2 
        delete(p);
        delete(q); 
    end
    p = plot3(ax3, x(1,:), x(2,:), x(3,:), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r', 'DisplayName', 'Trajectory');
    q = plot3(ax3, x(1,end), x(2,end), x(3,end), 'o', 'MarkerFaceColor', 'b', 'DisplayName', '$x_{\infty}$');
    pause(0.02);
end
legend(ax3, 'Interpreter', 'Latex');
hold off;
title(ax3, 'Animation Stable Invariant 3D');
%% 1c) Family of Solutions for Stable Node 3D

delete(p);
delete(q);
delete(initial_point);
delete(legend);

% in the 3d case we take some points on a sphere, therefore we need 2
% parameters theta and phi
phi = linspace(0, pi, 5);
theta = linspace(0, 2*pi, 10);

hold on;
delete(legend);
for i=1:length(theta)
    for j=1:length(phi)
        x0 = 7*[sin(phi(j))*cos(theta(i)); sin(phi(j))*sin(theta(i)); cos(phi(j))];
        x = solveLinearSystem(A,t,x0);
       
        if abs(phi(j)-pi/2)<1e-3 
            c = 'r'; 
        else 
            c='k'; 
        end
        plot3(ax3, x0(1,1), x0(2,1), x0(3,1),  Marker='o', Color=c, DisplayName='$x_0$', HandleVisibility='off'); %this is the initial point
        plot3(ax3, x(1,:), x(2,:), x(3,:), c, 'linewidth', 2, LineStyle='--', HandleVisibility='off');
    end
end
title(ax3, '3D stable node');
hold off;

%% 2) 3D Case - 2 stable and 1 unstable
clear; clc; close all;

A = [-5, 0, 0;
    0, -1, 0;
    0, 0, 1
    ];

[V, D] = eig(A);

disp('Eigenvalues of A:');
disp(D);

L = 8;
l = linspace(-L, L, 10).*sqrt(2); %range for the eigenspaces
v1 = ( V(:,1)./norm(V(:,1)) ) * l; % normalized eigenvector times l
v2 = ( V(:,2)./norm(V(:,2)) ) * l;
v3 = ( V(:,3)./norm(V(:,3)) ) * l;

figure;
ax4 = gca;
plot3(ax4, v1(1,:), v1(2,:), v1(3,:), ...
        v2(1,:), v2(2,:), v2(3,:), ...
        v3(1,:), v3(2,:), v3(3,:), 'linewidth', 2);
leg = legend(ax4, 'v1', 'v2', 'v3');
xlabel('x1');
ylabel('x2');
zlabel('x3');
grid(ax4);
axis(1.1 * [-L L -L L -L L]); 
legend(ax4,'v1','v2','v3 - Unstable', 'Interpreter', 'latex');

patch([v1(1,1) v2(1,1) v1(1,end) v2(1,end)],...
     [v1(2,1) v2(2,1) v1(2,end) v2(2,end)],...
     [v1(3,1) v2(3,1) v1(3,end) v2(3,end)],'b','FaceAlpha',0.2, 'HandleVisibility', 'off');

%% 2a) Animate the plot for 2 stable and 1 unstable eigenvectors

x0 =  5*V(:,1)+ -6*V(:,2) + 0.1*V(:,3);
eps = 1/50;
Tf = abs(log(norm(x0) / eps) / -max(D(D~= 0))); % where D are eigenvalues
t = 0:0.05:Tf;  % time span

hold on;
% Plot the initial point
delete(legend);
initial_point = plot3(ax4, x0(1), x0(2), x0(3), 'Marker','o', 'Color', 'black', 'DisplayName', '$x_0$');
p = [];
q = [];
for i = 2:length(t)
    x = solveLinearSystem(A, t(1:i), x0);
    if i > 2 
        delete(p);
        delete(q); 
    end
    p = plot3(ax4, x(1,:), x(2,:), x(3,:), 'LineWidth', 2, 'LineStyle', '--', 'Color', 'r', 'DisplayName', 'Trajectory');
    q = plot3(ax4, x(1,end), x(2,end), x(3,end), 'o', 'MarkerFaceColor', 'b', 'DisplayName', '$x_{\infty}$');
    pause(0.02);
end
legend(ax4, 'Interpreter', 'Latex');
hold off;
title(ax4, 'Animation 2 stable and 1 unstable eigenvectors 3D');

%% 2b) Family of Solutions for 2 stable and 1 unstable 3D

delete(p);
delete(q);
delete(initial_point);
delete(legend);

% in the 3d case we take some points on a sphere, therefore we need 2
% parameters theta and phi
phi = linspace(0, pi, 5);
theta = linspace(0, 2*pi, 10);

hold on;
delete(legend);
for i=1:length(theta)
    for j=1:length(phi)
        x0 = 7*[sin(phi(j))*cos(theta(i)); sin(phi(j))*sin(theta(i)); cos(phi(j))];
        x = solveLinearSystem(A,t,x0);
       
        if abs(phi(j)-pi/2)<1e-3 
            c = 'r'; 
        else 
            c='k'; 
        end
        plot3(ax4, x0(1,1), x0(2,1), x0(3,1),  Marker='o', Color=c, DisplayName='$x_0$', HandleVisibility='off'); %this is the initial point
        plot3(ax4, x(1,:), x(2,:), x(3,:), c, 'linewidth', 2, LineStyle='--', HandleVisibility='off');
    end
end
title(ax4, '3D unstable node');
hold off;

%% 3) 3D with complex eigenvalues and 1 real eigenvalue
clc; clear; close all;

lambda = -1;
alpha = - 0.4;
omega = 4;

% A = [lambda, 0, 0;
%     0, alpha, omega;
%     0, - omega, alpha
%     ];

% construct a matrix A having lambda_1,2 = sigma +- i omega and lambda_3 =
% lambda
V = [1 1 0;
    1 -1 1;
    -1/8 1/8 8];

A = V*[alpha omega 0;-omega alpha 0;0 0 lambda]*inv(V);

[V_complex, D] = eig(A);

disp('Matrix A');
disp(A);

disp('Eigenvalues of A:')
disp(D);

L = 8;
l = linspace(-L, L, 10).*sqrt(2); %range for the eigenspaces
v1 = ( V(:,1)./norm(V(:,1)) ) * l; % normalized eigenvector times l
v2 = ( V(:,2)./norm(V(:,2)) ) * l;
v3 = ( V(:,3)./norm(V(:,3)) ) * l;

figure;
ax5 = gca;
plot3(ax5, v1(1,:), v1(2,:), v1(3,:), ...
        v2(1,:), v2(2,:), v2(3,:), ...
        v3(1,:), v3(2,:), v3(3,:), 'linewidth', 2);
leg = legend(ax5, 'v1', 'v2', 'v3');
xlabel('x1');
ylabel('x2');
zlabel('x3');
grid(ax5);
axis(1.1 * [-L L -L L -L L]); 
legend(ax5,'v1','v2','v3', 'Interpreter', 'latex');

patch([v1(1,1) v2(1,1) v1(1,end) v2(1,end)], ...
     [v1(2,1) v2(2,1) v1(2,end) v2(2,end)], ...
     [v1(3,1) v2(3,1) v1(3,end) v2(3,end)],'b','FaceAlpha',0.2, 'HandleVisibility', 'off');
%% 3b) Animate the Solution
x0 =  5*V(:,1)./norm(V(:,1)) + 5*V(:,2)/norm(V(:,2)) + 6*V(:,3)/norm(V(:,3));

eps = 1/50;
Tf = 10;
t = 0:0.02:Tf;


p = [];
q = [];
hold on;
% Plot the initial point
delete(legend);
initial_point = plot3(ax5, x0(1), x0(2), x0(3), 'ko', 'MarkerFaceColor', 'k', DisplayName='$x_0$');
for i = 2:length(t)
    x = solveLinearSystem(A, t(1:i), x0);
    if i > 2 
        delete(p); 
        delete(q); 
    end
    p = plot3(ax5, x(1,:), x(2,:), x(3, :), 'LineWidth', 1.5,LineStyle='--', DisplayName='Trajectory', Color='red');
    q = plot3(ax5, x(1,end), x(2,end), x(3, end), 'bo', 'MarkerFaceColor', 'b', DisplayName='$x_{\infty}$');
    pause(0.02);
end

if alpha<0
    title(ax5, 'Stable focus');
else
    title(ax5, 'Unstable focus');
end
legend(ax5, 'Interpreter', 'Latex');
hold off;
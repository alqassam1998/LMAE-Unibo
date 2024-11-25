clc; clear all; close all;
% Analysis of Discrete Time Linear Systems
% x_(k+1) = A x_k

% Step 1: Specify the desired eigenvalues
lambda = [-1.5, -0.5]; % Example eigenvalues
n_dim = length(lambda);

% Step 2: Create a diagonal matrix with these eigenvalues
Lambda = diag(lambda);

% Step 3: Generate a random invertible matrix V
V = rand(n_dim) * 1i; % Generate a random 2x2 matrix

% Step 4: Construct the matrix A
A = V * Lambda * inv(V);

% A = [1 1; 0 1];

% Display the matrix A
disp('Matrix A:');
disp(A);

%%

% Verify that the eigenvalues of A match the specified values
eigenvalues = eig(A);
disp('Eigenvalues of A:');
disp(eigenvalues);

% Initial state
x0 = [1; 2];

% Number of time steps
numSteps = 30;

% Initialize state storage
x = zeros(2, numSteps);
x(:, 1) = x0;

% Simulate the system
for k = 2:numSteps
    x(:, k) = A * x(:, k-1);
end

figure;
ax1 = gca;
hold(ax1, 'all');
theta = linspace(0, 2*pi, 100); 
x_unit_circle = cos(theta); % x-coordinates for unit circle
y_unit_circle = sin(theta); % y-coordinates for unit circle

% Plot the unit circle
plot(ax1, x_unit_circle, y_unit_circle, 'b-', 'LineWidth', 1.5); 
% Plot the eigenvalues
plot(ax1, real(lambda), imag(lambda), 'rx', 'MarkerSize', 15, 'DisplayName', 'Eigenvalues'); % Red dots for eigenvalues
xlabel('Real Part');
ylabel('Imaginary Part');
title('Eigenvalues on the Unit Circle in the Gaussian Plane');
axis equal; 
xline(0, '--k'); 
yline(0, '--k'); 
legend('Unit Circle', 'Eigenvalues');
grid on;
xlim([-2 2]); 
ylim([-2 2]); 
hold off;

% Plot the results as discrete points
figure;
ax2 = gca;
hold(ax2, 'all');
for i=1:n_dim
    stem(ax2, 1:numSteps, x(i, :), 'filled', 'DisplayName', sprintf('x_{%d}', i), 'Marker', 'o');
end
xlabel('Time step (k)');
ylabel('State value');
title('Modal Response of the Discrete-Time System');
legend;
grid on

%%
clear all; close all; clc;
% Zero Eigenvalue

A = [0 1 0; 0 0 1; 0 0 0];
B = [0; 0; 1];
u = 0;

n_dim = length(A(:, 1));
numSteps = 30;

% Initialize state storage
x = zeros(n_dim, numSteps);
x0 = randi(10, [n_dim, 1]);
x(:, 1) = x0;

% Simulate the system
for k = 2:numSteps
    x(:, k) = A * x(:, k-1) + B * u;
end

eigenvalues = eig(A);
disp('Eigenvalues of A:');
disp(eigenvalues);

figure;
ax3 = gca;
hold(ax3, 'all');
stem(1:numSteps, x(1, :), 'filled', 'DisplayName', 'x_1', 'Marker', 'o');
% stem(1:numSteps, x(2, :), 'filled', 'DisplayName', 'x_2', 'Marker', 'o');
% stem(1:numSteps, x(3, :), 'filled', 'DisplayName', 'x_3', 'Marker', 'o');
xlabel('Time step (k)');
ylabel('State value');
title('Zero Eigenvalue Example - dynamics of a buffer');
legend;
grid on

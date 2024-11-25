%% SYSTEM THEORY AND ADVANCED CONTROL - Academic Year 2024/2025
% EXERCISE SESSION V
% Kalman Reachability Form, Control Canonical Form and Pole Placement

close all; clear all; clc;

% The system is a Discrete Time one

A = [0 0 0 0; % n = 4
    0 0 0 0;
    0 1 0 0;
    1 0 0 1
    ];

B = [0 0;
    0 1;
    0 0;
    1 0
    ];

R = ctrb(A, B)
rank(R)
%%

T_1 = [R(:, 1), R(:, 4), R(:, 2), [1; 0; 0; 0]]

%%

A_kalman = pinv(T_1) * A * T_1

B_kalman = T_1 * B

%%
% Simulation of the Free Response
x0 = rand(size(A, 1), 1);
n_steps = 20;
xx = zeros(size(x0, 1), n_steps);
xx(:, 1) = x0;
 for i=2:(n_steps + 1)
     xx(:, i) = A_kalman * xx(:, i - 1) + B_kalman * [0; 0];
 end



% Does the eigenvalues of the Kalman decomposition coincide with the
% original matrix?
disp('Eigenvalues of the Kalman decomposition')
disp(eig(A_kalman))

disp('Eigenvalues of the Original Matrix')
disp(eig(A))

% Is this system Stabilizable? The NR submatrix is Schur, it is.
A_r = A_kalman(1:3, 1:3);
B_r = B_kalman(1:3, 1:2);

disp('Reachable Submatrix Ar')
disp(A_r)

%%
% Desired Eigenvalues
desired_eig = [-0.1, -0.2, 0];
K_r = place(A_r, - B_r, desired_eig); 
% using Place the standard notation is A - BK, -Br is done to match ours. 
% Moreover, the multiplicity of the eigenvalues that you can set is linked
% to the number of inputs (in this case, you can not put 3 eigenvalues at
% zero because we have only two inputs).

n_nr = size(A, 1) - size(A_r, 1);
K_kalman = [K_r, zeros(size(B, 2), n_nr)];

K = K_kalman * pinv(T_1);

% Back to the original System
A_feedback = (A + B * K);

disp('Desired Eigenvalues')
disp(desired_eig)

disp('Eigenvalues of A + BK')
eig(A_feedback)

x0 = rand(size(A, 1), 1);
n_steps = 10;
xx = zeros(size(x0, 1), n_steps);
xx_aft_feed = zeros(size(x0, 1), n_steps);
xx(:, 1) = x0;
xx_aft_feed(:, 1) = x0;

% Simulation Loop
 for i=2:(n_steps + 1)
     xx(:, i) = A_kalman * xx(:, i - 1) + B_kalman * [0; 0];
     xx_aft_feed(:, i) = A_feedback * xx_aft_feed(:, i - 1);
 end

figure;
ax = gca;
hold(ax, 'on');
title(ax, 'Free Response Kalman Coordinate');
stairs(ax, 0:n_steps', xx(1, :)', LineWidth=1.5, DisplayName='$x_1$');
stairs(ax, 0:n_steps', xx(2, :)', LineWidth=1.5, DisplayName='$x_2$');
stairs(ax, 0:n_steps', xx(3, :)', LineWidth=1.5, DisplayName='$x_3$');
stairs(ax, 0:n_steps', xx(4, :)', LineWidth=1.5, DisplayName='$x_4$');
legend(ax, interpreter="latex")
xlabel('Steps')
ylabel('System States')
grid(ax)

figure;
ax2 = gca;
hold(ax2, 'on');
title(ax2, 'Behavior After Feedback');
stairs(ax2, 0:n_steps', xx_aft_feed(1, :)', LineWidth=1.5, DisplayName='$x_1$');
stairs(ax2, 0:n_steps', xx_aft_feed(2, :)', LineWidth=1.5, DisplayName='$x_2$');
stairs(ax2, 0:n_steps', xx_aft_feed(3, :)', LineWidth=1.5, DisplayName='$x_3$');
stairs(ax2, 0:n_steps', xx_aft_feed(4, :)', LineWidth=1.5, DisplayName='$x_4$');
legend(ax2, interpreter="latex")
xlabel('Steps')
ylabel('System States')
grid(ax2)

%% Controllability Canonical Form
% Consider the CT system
clc; clear all; close all;
A = [1 2 0;
    0 0 1;
    0 1 0
    ];

B = [1; 0; 1];

R = ctrb(A, B);
rank(R)
det(R)
%%
% Define symbolic variables
syms lambda
A = sym(A);        

% Compute the characteristic polynomial
I = eye(size(A));              
char_poly = det(A - lambda*I);
char_poly2 = det(lambda*I - A);

% Display the characteristic polynomial
disp('The characteristic polynomial is:')
disp(char_poly)
disp(char_poly2)

coeffs_char_poly2 = coeffs(char_poly2, lambda, 'All');
degree_charact = length(coeffs_char_poly2) - 1;

coeffs_char_poly2 = coeffs_char_poly2 / coeffs_char_poly2(1); % ensure that lambda^n is monic

alphas = fliplr(- coeffs_char_poly2(2:end));

%%
A_con = [zeros(degree_charact - 1, 1), eye(degree_charact - 1); alphas]
B_con = [0; 0; 1]

disp('System Matrices in Controllability Form');
disp(A_con), disp(B_con);

% It is possible to compute it also using the relation T = R_c R-1
A = [1 2 0;
    0 0 1;
    0 1 0
    ];

R_c = ctrb(A_con, B_con);
R_1 = pinv(R);

Tc = R_c * R_1;

A_con2 = Tc * A * pinv(Tc);
B_con2 = Tc * B;

disp('Controllability Matrices Using T');
disp(A_con2), disp(B_con2);
%%
% Let's compare the result obtained by hand and the place command
% We want poles of the original system in [-2, -2, -5]
desired_poles = [-2, -2, -5];
K_c = [-19, -25, -10]; % Obtained by solving gamma(A + BK) = gamma(A_c + B_cK_c)
K = K_c * Tc;

A_feedback_hand = (A + B*K);
disp('Eigenvalues of A + BK computed by hand')
disp(eig(A_feedback_hand))

% Using place
% desired_poles = [-2, -5, -10]; % DO NOT RELY ON PLACE
% K_place_c = place(A_con, - B_con, desired_poles); %%% !!!!! Does not mean that you could not do that. Place is not able to.
% K_place = K_place_c * Tc;
% A_feedb_place = (A + B*K_place);
% disp('Eigenvalues of A + BK using place')
% disp(eig(A_feedb_place))


% Simulation
model_feedback = @(t, x) [A_feedback_hand * x];
x0 = 10 * rand(size(A, 1), 1);
dt = 1e-4;
t_end = 5;
[time, xx] = ode45(model_feedback, 0:dt:t_end, x0);

figure;
ax1 = gca;
hold(ax1, "on")
plot(ax1, time, xx, LineWidth=1.5)
xlabel('Time (s)', Interpreter='latex')
ylabel('System states', Interpreter='latex')
legend(ax1, '$x_1$', '$x_2$', '$x_3$', interpreter='latex')
grid(ax1)
hold(ax1, "off");

%% Exercise
% State Space Model
A = [-3.0   0.0  0.0;
      0.0   0.0  1.0;
      0.0  -2.0  3.0];
  
B = [0.0;
     0.0;
     1.0];

C = [1.0, 0.0, 0.0;
     0.0, 1.0, 0.0];

% Check for Stability
% --> Your Code Here

% Check for Controllability
% --> Your Code Here

% PBH Test for Controlllability
% --> Your Code Here

% Compute the Kalman Controllability Decomposition
% --> Your Code Here

% Show the Reachable/Not Reachable Submatrix
% --> Your Code Here

% Check Stability of the Reachable Part
% --> Your Code Here

% Design K in Kalman Coordinates for the Reachable Part Only
poles = [-5, -6];

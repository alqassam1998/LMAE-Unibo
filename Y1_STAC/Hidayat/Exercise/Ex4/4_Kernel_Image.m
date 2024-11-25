function [kernel_basis, image_basis, image_dimension] = compute_kernel_and_image(A)
    % Input:
    %   A - Input matrix
    % Output:
    %   kernel_basis   - Basis for the kernel (null space) of A, with unit length vectors
    %   image_basis    - Basis for the image (column space) of A, with unit length vectors
    %   image_dimension - Dimension of the image space (rank of A)

    % Calculate the kernel (null space) of A
    kernel_basis = null(A);  % Basis for the null space of A (not orthonormal)

    % Normalize each vector in the kernel basis to unit length
    for i = 1:size(kernel_basis, 2)
        kernel_basis(:, i) = kernel_basis(:, i) / norm(kernel_basis(:, i));
    end

    % Find independent columns for the image (column space)
    [~, pivot_columns] = rref(A);  % Identify pivot columns in reduced row echelon form
    image_basis = A(:, pivot_columns);  % Select columns corresponding to pivots

    % Normalize each vector in the image basis to unit length
    for i = 1:size(image_basis, 2)
        image_basis(:, i) = image_basis(:, i) / norm(image_basis(:, i));
    end

    % Calculate the dimension of the image space (rank of A)
    image_dimension = size(image_basis, 2);  % Number of pivot columns is the rank
end


% Define samples matrices
C = 1e-1;
L1 = 1;
L2 = 2;

A = [0, 0, - 1/L1;
    0, 0, -1/L2;
    1/C, 1/C, 0];

B = [1/L1; 1/L2; 0
    ];

C = [1, 1, 0];

O = obsv(A, C);
rank(O);

% Compute orthonormal bases for kernel and image
[~, image_basis, image_dimension] = compute_kernel_and_image(O');

% Display the results

disp('Orthonormal basis for the image (column space) of O^T:');
disp(image_basis);

disp(image_dimension);

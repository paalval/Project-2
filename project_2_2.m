clear all

rho_0 = 0; rho_n = 10; n = 1e2; % Start and endpoints + Number of mesh points
h = (rho_n - rho_0)/n; % Step length

for i=1:n
    rho(i) = rho_0 + i*h; % rho array
end

V = rho.^2; % Harmonic oscillator potential

e = (-1/h^2)*ones(1, n-1);                % Off-diagonal elements
omega_r = [0.01 0.5 1 5];                 % reflects strength of potential

fprintf('Two interacting electrons in harmonic oscillator potential eigenvalue problem\n')

legendCounter = 1;

for j = 1:numel(omega_r)

    fprintf('..........\n')

%     d = (2/h^2) + (omega_r(j)^2).*V + 1./rho; % diagonal entries - two interacting electrons
    d = (2/h^2) + (omega_r(j)^2).*V;          % diagonal entries - two non-interacting electrons

    % Tridiagonal matrix A - represents the Hamiltonian for the harmonic
    % oscillator
    A = gallery('tridiag', e, d, e);
    A = full(A); % Convert from sparse matrix to full matrix.
    B = A; % used for Jacobi; A for reference.

    % Calculate eigenvalues of A - these are the permitted energies in the
    % oscillator

    fprintf('for %d-by-%d symmetric tridiagonal matrix A\n', n, n)
    fprintf('The potential parameter is %0.2f\n', omega_r(j))

    fprintf('..........\n')

    fprintf('Jacobi method\n')

    % Use the Jacobi algorithm

    tolerance = 1E-8; % Tolerance
    i = 1; max_i = 1E8;
    [k, l] = offdiag(B);
    maxnondiag = abs(B(k, l));
    R = eye(n);
    T = R;
    while (maxnondiag > tolerance) && (i <= max_i)
        % Find max-valued non-diagonal element in B
        [k, l] = offdiag(B);
        maxnondiag = abs(B(k, l));

        % Rotate matrix B
        [B, T] = Jacobi(B, T, k, l, n);

        i = i + 1;
    end

    fprintf('Need %d transformations\n', i)

    fprintf('..........\n')

    % Test if orthogonality of the eigenvectors is conserved
    TT = (transpose(T)*T > tolerance); kroenicker_delta = eye(n);
    if isequal(TT, kroenicker_delta) == 1
        fprintf('Orthogonality of eigenvectors conserved\n');
    else
        fprintf('Orthogonality of eigenvectors not conserved\n')
        return % Stops the program
    end

    lambda_Jacobi = diag(B);
    [lowest1, index1] = min(lambda_Jacobi);  % Lowest eigval with position
    ground1 = T(:, index1); % Eigvec corresponding to lowest eigval

    % Plot the probability distribution of the ground state
    figure(1)
    plot(rho, ground1.^2')
    hold on

    fprintf('..........\n')

    % Use MATLAB's own function for finding eigenvalues

    fprintf('MATLAB eigenvalue function\n')

    [eigvec, lambda_MATLAB] = eig(A);
    [lowest2, index2] = min(diag(lambda_MATLAB));  % Lowest eigval with position
    ground2 = eigvec(:, index2); % Eigvec corresponding to lowest eigval

    % Plot the probability distribution of the ground state
    figure(2)
    plot(rho, ground2.^2)
    hold on
    
    % Legend strings
    ar{legendCounter} = sprintf('$\\omega_r = %0.2f$', omega_r(j));
    legendCounter = legendCounter + 1;
end

figure(1)
xlabel('Relative distance')
ylabel('Probability distribution of the ground state')
title('$f(\rho$) - Jacobi', 'interpreter', 'latex')
grid on
legend(ar, 'interpreter', 'latex')

figure(2)
xlabel('Relative distance')
ylabel('Probability distribution of the ground state')
title('$f(\rho$) - MATLAB', 'interpreter', 'latex')
grid on
legend(ar, 'interpreter', 'latex')
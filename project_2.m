clear all

rho_0 = 0; rho_n = 10; n = 1e2; % Start and endpoints + Number of mesh points
h = (rho_n - rho_0)/n; % Step length

for i=1:n
    rho(i) = rho_0 + i*h; % rho array
end

V = rho.^2; % Harmonic oscillator potential

e = (-1/h^2)*ones(1, n-1); % Off-diagonal elements
d = (2/h^2) + V;           % diagonal entries - single electron

% Tridiagonal matrix A - represents the Hamiltonian for the harmonic
% oscillator
A = gallery('tridiag', e, d, e);
A = full(A); % Convert from sparse matrix to full matrix.
B = A; % used for Jacobi; A for reference.

% Calculate eigenvalues of A - these are the permitted energies in the
% oscillator

fprintf('Single electron in harmonic oscillator potential eigenvalue problem\n')
fprintf('for %d-by-%d symmetric tridiagonal matrix A\n', n, n)

fprintf('..........\n')

% Use the Jacobi algorithm and time the process

tic

tolerance = 1E-8; % Tolerance
i = 1; max_i = 1E8;
[k, l] = offdiag(B);
maxnondiag = abs(B(k, l));
R = eye(n);
% R = eye(3);
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

fprintf('..........\n')

fprintf('Jacobi method\n')

toc

lambda_Jacobi = diag(B);

for k = 1:3
    text = 'lambda(%d) = %0.6f\n';
    fprintf(text, k, lambda_Jacobi(k))
end

fprintf('..........\n')

% Use MATLAB's own function for finding eigenvalues and time the process
tic

lambda_MATLAB = eig(A);

fprintf('MATLAB eigenvalue function\n')

toc

for k = 1:3
    text = 'lambda(%d) = %0.6f\n';
    fprintf(text, k, lambda_MATLAB(k))
end
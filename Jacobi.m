function [M, T] = Jacobi(A, R, k, l, n)
    % Jacobi algorithm for square n-by-n matrix A. Takes in a square matrix
    % A, indecies k and l, and dimensionality of the matrix n. Performs a
    % Jacobi rotation on A and returns the rotated matrix M. Matrix R is a
    % matrix containing all the eigenvectors as the column (used to check
    % if orthogonality is conserved).
    if A(k, l) ~= 0
        tau = (A(l, l) - A(k, k))/(2 * A(k, l));
        
        if tau >= 0
            t = -tau + sqrt(1 + tau*tau);
        else
            t = -tau - sqrt(1 + tau*tau);
        end
        
        c = 1/sqrt(1 + t*t);
        s = c*t;
    else
        c = 1;
        s = 0;
    end
    
    a_kk = A(k ,k);
    a_ll = A(l, l);
    A(k, k) = c*c*a_kk - 2*c*s*A(k, l) + s*s*a_ll;
    A(l, l) = s*s*a_kk + 2*c*s*A(k, l) + c*c*a_ll;
    A(k, l) = 0;    % Hard-coding off-diag elements by hand
    A(l, k) = 0;    % Ditto!
    
    for i = 1:n
        if (i ~= k) && (i ~= l)
            a_ik = A(i, k);
            a_il = A(i, l);
            A(i, k) = c*a_ik - s*a_il;
            A(k, i) = A(i, k);
            A(i, l) = c*a_il + s*a_il;
            A(l, i) = A(i, l);
        end
        % And finally the new eigenvectors
        r_ik = R(i,k);
        r_il = R(i,l);

        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    end
    M = A;
    T = R;
end


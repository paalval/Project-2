function [k, l] = offdiag(A)
    % Takes in a symmetric (real) matrix A and returns the indecies for
    % the maximum value on the off-diagonal in matrix A.
    % Example:
    %   A = [1 4 3; 8 9 2];
    %   [k, l] = offdiag(A);
    % Will give row k = 2 and column l = 1. A(2, 1) now gives 8, which is
    % the largest off-diagonal element in the example matrix, 
    
    % Determines the dimensions of matrix A.
    [n, m] = size(A);
    
    % Set a starting point for each iteration over j to consider
    max_val = 0;
    
    % Determine indecies of maximum off-diagonal element
    for i = 1:n
        for j = (i + 1):m
            a_ij = abs(A(i, j))^2;
            if a_ij > max_val;
                max_val = a_ij; % Set new maximum
                k = i;          % Row index
                l = j;          % Column index
            end
        end
    end
end


function correlated_matrix = Shadowing_matrix(rows, sigma_s, lambda)
    % Generates a correlated matrix based on input parameters and extracts values
    % corresponding to coordinates in SU.
    % 
    % Parameters:
    % rows: Number of rows (and columns) in the correlation matrix.
    % sigma_s: Standard deviation of the shadowing in dB.
    % SU: Coordinates of m points randomly distributed in a circular area.
    %
    % Returns:
    % correlated_values: Extracted values from the correlated matrix.

    % Generate uncorrelated Gaussian (normal) random matrix
    uncorrelated_matrix = randn(rows);

    % Form the distance matrix
    [X, Y] = meshgrid(1:rows, 1:rows);
    distances = sqrt((X - X').^2 + (Y - Y').^2);

    % Compute the correlation matrix based on the negative-exponential model
    correlation_matrix = exp(-distances / lambda);

    % Cholesky decomposition to get the lower triangular matrix
    L = chol(correlation_matrix, 'lower');

    % Multiply the uncorrelated matrix by the Cholesky factor
    correlated_matrix = L * uncorrelated_matrix * L';

    % Rescale the matrix to have unit variance
    correlated_matrix = sigma_s * correlated_matrix;%/std(correlated_matrix(:));
    % correlated_matrix = flipud(correlated_matrix');

   end

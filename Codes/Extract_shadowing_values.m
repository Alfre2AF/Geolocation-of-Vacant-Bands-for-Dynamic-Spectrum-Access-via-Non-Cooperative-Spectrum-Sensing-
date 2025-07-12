function m_shadowing_values = Extract_shadowing_values(SU, correlated_matrix, r)
    
    rows = size(correlated_matrix,1);

    % Normalize the SU coordinates from [-r, r] to [1, rows]
    SU = (SU + r) * (rows - 1) / (2 * r) + 1;

    % Round the coordinates to the nearest integer for indexing
    SU = round(SU);

    % Ensure the indices are within the bounds of the matrix
    SU(SU < 1) = 1;
    SU(SU > rows) = rows;

    % Extract the values from the correlated_matrix at the specified coordinates
    m_shadowing_values = correlated_matrix(sub2ind([rows, rows], SU(:, 1), SU(:, 2)));
end

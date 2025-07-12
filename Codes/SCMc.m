function cov_matrix = SCM(X, M)
    % Calcular la media de X
    mu = mean(X);
    % Centramos el vector X restando la media
    X_centered = X - mu;
    
    % Longitud de X_centered
    n = length(X_centered);
    
    % Verificar que M divida n
    if mod(n, M) ~= 0
        error('M must divide n without a remainder.');
    end
    
    % Número de sub-vectores
    num_sub_vectors = n / M;
    
    % Inicializar matriz para almacenar sub-vectores
    sub_vectors = zeros(M, num_sub_vectors);
    
    % Crear sub-vectores no solapados
    for i = 1:num_sub_vectors
        start_idx = (i-1)*M + 1;
        sub_vectors(:, i) = X_centered(start_idx:start_idx+M-1);
    end
    
    % Inicializar matriz de covarianza
    cov_matrix = zeros(M, M);
    
    % Calcular la covarianza promediando los productos exteriores
    for i = 1:num_sub_vectors
        deviation = sub_vectors(:, i); % Ya centrado, media = 0
        cov_matrix = cov_matrix + (deviation * deviation');
    end
    
    % Normalizar
    cov_matrix = cov_matrix / num_sub_vectors;
end

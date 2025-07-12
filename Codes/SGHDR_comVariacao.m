
%% Abordagem do Problema de Geolocalização de Buracos de Espectro para Acesso Dinâmico ao Espectro
%% Variacion de parametros
%% Alfredo Jesus Arbolaez Fundora.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables; clc;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parâmetros do sistema

N1 =1600;              % Número de nós SSIoT na rede base.
m = 1;
G = 4;
M = 10;
L = 2;                 % Comprimento da área de cobertura quadrada centrada em (0,0).
alpha = 0.3;           % Fração de SSIoTs fixos (alpha*N deve ser inteiro).
n = 2400;              % Número de amostras por SU. Quando m=1, defina n igual ao número total de amostras no CSS.
SNR = -5;              % Relação sinal-ruído média sobre todos os SSIoTs, em dB.
runs = 5;          % Número de rodadas de sensoriamento para computar as CDFs empíricas.
eta = 2.5;             % Expoente de perda de percurso.
d0 = 0.001*L;          % Distância de referência para o cálculo da perda de percurso, em metros.
P_tx = 5;              % Potência de transmissão do PU, em watts.
xPU = 3*L;             % Coordenada x=y do transmissor PU, em metros. Quando distante dos SSIoTs => área de proteção grande para os PUs.
Ns = 3*M;              % Número de amostras por símbolo QPSK do PU.
rho = 0.5;             % Fração das variações de potência do ruído em torno da média.
meanK = 1.88;          % Média do fator Rice (em dB) para K variável ao longo das rodadas e dos SSIoTs.
sdK = 4.13;            % Desvio padrão (em dB) de K ao longo das rodadas e dos SSIoTs.
randK = 1;             % Se randK = 1, K é aleatório; se randK = 0, K = meanK.
PUsignal = 1;          % Sinal do PU: "0" = Gaussiano niid; "1" = QPSK niid (Ns>1) ou iid (Ns=1).
Cc = 0.95;              % Coeficiente de correlação entre amostras vizinhas do sinal Gaussiano niid.
sigma_s = 7;           % Desvio padrão do sombreamento, em dB.
rows = 50;             % Número de linhas da matriz completa de sombreamento. Veja a função Shadowing_matrix.m.
Lambda = 0.8*rows;     % Comprimento de correlação (Lambda > 0 é uma fração de 'rows').
Npt = 100;             % Número de pontos nas ROC.
SNR_th = SNR-7;          % Limiar de SNR (em dB) que define as áreas sombreadas.
P_rx_th = -55;
Pfa = 0.25;
Factor = 4.25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here, enable only the set of parameter values to be varied:

%Parameter = [500, 800, 1000, 1200, 1600, 2000];  Flag=1; % For varying N1.
% Parameter = [5, 10, 15, 20, 25, 30];          Flag=2; % For varying M.
%Parameter = [-11,-8,-5,-2,0,2];             Flag=3; % For varying SNR.
%Parameter = [1000, 1500, 2000, 2400, 3000, 4000];  Flag=4;
Parameter = [3.25, 3.5, 3.75, 4,4.25, 4.5 ];  Flag=5;

for par = 1 : size(Parameter,2) % VARIABLE PARAMETER LOOP
 
if Flag==1; N1 = Parameter(par); end

if Flag==2; M = Parameter(par); end


if Flag==3; SNR = Parameter(par); end

if Flag==4; n = Parameter(par); end

if Flag==5; Factor = Parameter(par); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Valor Esperado de P(d) via integração numérica (instabilidade numérica pode ocorrer se o PU estiver dentro da área de cobertura)
integrand = @(x, y) P_tx.*(sqrt((x - xPU).^2 + (y - xPU).^2)/d0).^-eta;
E_Pd = integral2(integrand, -L/2, L/2, -L/2, L/2) / (L^2);

%% Valor Esperado de P_rx
E_10S10 = exp(((sigma_s^2)*log(10)^2)/200); % Valor esperado de 10^(S/10)
E_Prx = E_Pd * E_10S10; % Valor esperado de P_rx (com sombreamento)

%% Variância média do ruído, de acordo com o SNR
if rho > 0
    Sigma2_avg = (E_Prx*(log((1+rho)/(1-rho)))/(2*rho))/(10^(SNR/10)); % Sigma2_avg para 0 < rho < 1.
else
    Sigma2_avg = E_Prx/(10^(SNR/10)); % Sigma2_avg = 1 para rho = 0.
end

%% Pré-alocação de variáveis
Tpride_h0 = zeros(runs,N1); 
Tpride_h1 = zeros(runs,N1);
snr = zeros(runs,N1); % Pré-alocação do vetor para armazenar as medidas de SNR

Local_decisions_H1H0 = zeros(runs,N1); 
Local_decisions_H1H1 = zeros(runs,N1);
P_fa = zeros(Npt,N1); 
P_d = zeros(Npt,N1);
Frac_SHG_correct = zeros(1,runs);
actual_SHG = zeros(rows, rows, runs);

%% Distribuição aleatória de alpha*N nós fixos dentro da área de cobertura dos SSIoTs
x1 = rand(1, alpha*N1)*L-L/2; 
y1 = rand(1, alpha*N1)*L-L/2;
%% Distribuição aleatória de (1-alpha)*N nós móveis dentro da área de cobertura dos SSIoTs
x2 = [x1 , rand(1, N1-alpha*N1)*L-L/2]; 
y2 = [y1 , rand(1, N1-alpha*N1)*L-L/2];

%% Eventos de sensoriamento para todos os nós
for run = 1:runs
    fprintf('Ejecutando run %d de %d...\n', run, runs);
    %% Valores de sombreamento para todos os pontos (rows^2) na área
    Full_shadowing_matrix = Shadowing_matrix(rows, sigma_s, Lambda);

    %% Distâncias do PU para todos os pontos (rows^2) na área
    d_all = zeros(rows,rows); 
    for j = 0:rows-1
        for u = 0:rows-1
            d_all(j+1,u+1) = norm([((j-rows/2)/(rows/2))*L/2, ((u-rows/2)/(rows/2))*L/2] - [xPU, xPU]); % Distâncias de (x,y) até (xPU,xPU)
        end
    end

    %% Potência do sinal recebido para todos os pontos (rows^2) na área (usado para calcular os buracos de espectro reais)
    P_rxdBm_all = 10*log10(10^3 * P_tx*(d0./d_all).^eta) + Full_shadowing_matrix; 
    P_rx_all = 10^-3 * 10.^(P_rxdBm_all/10);
    SNR_all = 10*log10(P_rx_all/Sigma2_avg);
    actual_SHGsnr = SNR_all < SNR_th;
  actual_SHG(:,:,run) = P_rxdBm_all < P_rx_th;


    %% Loop para os nós
    for nodos = 1:N1
        %% Coordenadas dos SSIoTs para um nó
        SSIoT = [x2(nodos), y2(nodos)];

        %% Valores de sombreamento para os SSIoTs de um nó
        Shadowing_values = Extract_shadowing_values(SSIoT, Full_shadowing_matrix, L/2);

        %% Distâncias do PU aos SSIoTs de um nó
        d_pu = norm(SSIoT - [xPU, xPU]);  % Apenas um escalar

        %% Potências dos sinais recebidos pelos SSIoTs de um nó
        P_rxdBm = 10*log10(10^3 * P_tx*(d0./d_pu).^eta) + Shadowing_values; 
        P_rx = 10^-3 * 10.^(P_rxdBm/10);
        g = sqrt(P_rx/P_tx); % Ganhos do canal dependentes da distância

        %% Variâncias do ruído (m x 1) variáveis ao longo das rodadas de sensoriamento
        U = unifrnd(-1, 1, m, 1); % RV uniforme em [-1,1].
        Sigma2 = (1 + rho*U)*Sigma2_avg; % Variâncias do ruído nos SSIoTs

        %% Vetor do canal (m x 1):
        a = zeros(m,1);
        for row = 1:m
            if randK == 1
                K = 10^((randn(1,1)*sdK+meanK)/10);
            else
                K = 10^(meanK/10); % K fixo
            end
            a(row,1) = (normrnd(sqrt(K/(2*(K+1))), sqrt((1-K/(K+1))/2), 1,1) + 1j*normrnd(sqrt(K/(2*(K+1))), sqrt((1-K/(K+1))/2), 1,1));
        end
        h = a.*g; % Redimensionar o vetor do canal de acordo com os ganhos acima

        %% Matrizes de ruído Gaussiano (m x n):
        W0 = zeros(m,n); 
        W1 = zeros(m,n);
        for j = 1:m
            W0(j,:) = normrnd(0, sqrt(Sigma2(j)/2), 1, n) + 1j*normrnd(0, sqrt(Sigma2(j)/2), 1, n);
            W1(j,:) = normrnd(0, sqrt(Sigma2(j)/2), 1, n) + 1j*normrnd(0, sqrt(Sigma2(j)/2), 1, n);
        end

        %% Sinal do PU (n x 1):
        if PUsignal == 0
            S = niid_Gaussian(n, P_tx, Cc); % Veja a função niid_Gaussian.m
        elseif PUsignal == 1
            numSymbols = ceil(n / Ns);
            realPart = (randi([0, 1], 1, numSymbols) * 2 - 1) .* ones(Ns, numSymbols);
            imagPart = (randi([0, 1], 1, numSymbols) * 2 - 1) .* ones(Ns, numSymbols);
            S = realPart + 1j * imagPart;
            S = S(:)'; % Achatar para um vetor linha
            S = circshift(S, randi([0, Ns-1])); % Deslocamento circular aleatório
            S = S(1:n) * sqrt(P_tx / 2); % Ajustar a potência
        end

        %% SNR medido em cada rodada para um nó
        snr(run,nodos) = mean(sum(abs(h*S).^2,2)./sum(abs(W0).^2,2));

        %% Matrizes de sinal recebido sob H0 e H1 (m x n)
        X_h0 = W0; 
        X_h1 = h*S + W1;

        %% Matrizes de covariância amostral (SCM)
        R_h0 = SCMc(X_h0,M); 
        R_h1 = SCMc(X_h1,M); % Veja a função SCM.m

        %% Estatística de teste PRIDe
        x_h0 = R_h0(:); 
        x_h1 = R_h1(:); 
        m0 = mean(x_h0); 
        m1 = mean(x_h1);
        Tpride_h0(run,nodos) = sum(abs(x_h0))/sum(abs(x_h0 - m0));
        Tpride_h1(run,nodos) = sum(abs(x_h1))/sum(abs(x_h1 - m1));
    end % Fim dos nós



end % Fim das rodadas (sensoriamento)

%% CDFs empíricas da estatística PRIDe para todos os nós
for nodos = 1:N1
    u = 0; 
    Th0 = Tpride_h0(:,nodos); 
    Th1 = Tpride_h1(:,nodos);
    Min = mean(Th0) - 3*std(Th0); 
    Max = mean(Th1) + 1*std(Th1);
    for i = Min:(Max-Min)/(Npt-1):Max
        aux_h0 = 0; 
        aux_h1 = 0; 
        u = u+1;
        for ii = 1:runs
            if Th0(ii) < i
                aux_h0 = aux_h0+1;
            end
            if Th1(ii) < i
                aux_h1 = aux_h1+1;
            end
        end
        CDF_Tpride_H0(u,nodos) = aux_h0/runs;
        CDF_Tpride_H1(u,nodos) = aux_h1/runs;
    end
end

%% Limiares de decisão dos nós e Pd para o Pfa alvo
for nodos = 1:N1
    Z = sort(Tpride_h0(:,nodos));
   Gamma(nodos) = Z(round((1-Pfa)*runs)); % limiar de decisão
    %Gamma(nodos) = prctile(Tpride_h0(:, nodos), 100*(1 - Pfa));

    aux_h0 = 0; 
    aux_h1 = 0;
    for ii = 1:runs
        if Tpride_h1(ii,nodos) >= Gamma(nodos)
            aux_h1 = aux_h1 + 1;
        end
        if Tpride_h0(ii,nodos) >= Gamma(nodos)
            aux_h0 = aux_h0 + 1;
        end
    end
    Pdc(nodos) = aux_h1/runs; % Apenas para depuração (comparar com as ROCs dos nós para o Pfa local alvo)
    Pfac(nodos) = aux_h0/runs; % Apenas para depuração (comparar com as ROCs dos nós para o Pfa local alvo)
    % Decisões dos nós para H1|H0 com o limiar de decisão dado
    Local_decisions_H1H0(:,nodos) = Tpride_h0(:,nodos) >= Gamma(nodos);
    % Decisões dos nós para H1|H1 com o limiar de decisão dado
    Local_decisions_H1H1(:,nodos) = Tpride_h1(:,nodos) >= Gamma(nodos);

    
end

%% Pfa e Pd locais a partir das CDFs empíricas
for nodos = 1:N1
    P_fa(:,nodos) = 1 - CDF_Tpride_H0(:,nodos);
    P_d(:,nodos)  = 1 - CDF_Tpride_H1(:,nodos);
end

%% Avaliação global a partir das decisões individuais.
% Supõe-se que Local_decisions_H1H0 e Local_decisions_H1H1 são matrizes de dimensão (runs x N1)
global_Pfa = zeros(runs, 1);
global_Pd  = zeros(runs, 1);
for run = 1:runs
    % Para cada rodada, calcula-se a porcentagem de sensores que tomam uma decisão de "detecção"
    global_Pfa(run) = sum(Local_decisions_H1H0(run, :)) / N1;
    global_Pd(run)  = sum(Local_decisions_H1H1(run, :)) / N1;
end

mean_global_Pfa = mean(global_Pfa);
mean_global_Pd  = mean(global_Pd);

disp(['Global Pfa (média): ', num2str(mean_global_Pfa)]);
disp(['Global Pd (média): ', num2str(mean_global_Pd)]);





%% Matriz extendida de decisión usando Voronoi e interpolación (por run)
for run = 1:runs
     x = x2;  % Armazena as posições finais de todos os sensores
    y = y2;
    % --- 1. Obtener las coordenadas de los sensores y sus decisiones ---
    % Se asume que 'x' y 'y' contienen las posiciones de los sensores,
    % y que 'Local_decisions_H1H1' almacena las decisiones locales para H1.
    X = x(:);
    Y = y(:);
    % Aquí se utiliza la decisión de la última ronda; si runs>1, puede cambiarse
    decision_nodos = Local_decisions_H1H1(run, :); 
    
    % --- 2. Calcular el diagrama de Voronoi ---
    [V, C] = voronoin([X, Y]);
    
    % --- 3. Interpolación entre nodos que no detectaron señal ---
    % Se define un umbral de distancia
    d_thresh = Factor / sqrt(N1 / (L^2));
    % Inicializar vector para marcar sensores (blanco) que se encuentren próximos
    marcar_blanco = false(length(X), 1);
    % Índices de nodos que NO detectaron señal
    idx_nodetect = find(decision_nodos == 0);
    % Para cada par de nodos que no detectaron, si la distancia es menor al umbral,
    % se marca la región entre ellos para forzar la detección (interpolación)
    for i = 1:length(idx_nodetect)
        for j = i+1:length(idx_nodetect)
            ni = idx_nodetect(i);
            nj = idx_nodetect(j);
            dist = hypot(X(ni) - X(nj), Y(ni) - Y(nj));
            if dist <= d_thresh
                x_min = min(X(ni), X(nj)); 
                x_max = max(X(ni), X(nj));
                y_min = min(Y(ni), Y(nj)); 
                y_max = max(Y(ni), Y(nj));
                en_rect = (X >= x_min & X <= x_max & Y >= y_min & Y <= y_max);
                marcar_blanco(en_rect) = true;
            end
        end
    end
     nodos_aislados = false(length(X), 1);
    for i = 1:length(idx_nodetect)
        ni = idx_nodetect(i);
        distancias = sqrt((X(ni) - X(idx_nodetect)).^2 + (Y(ni) - Y(idx_nodetect)).^2);
        distancias(distancias == 0) = inf;
        if all(distancias > d_thresh)
            nodos_aislados(ni) = true;
        end
    end
    % --- 4. Crear la malla para la interpolación ---
    % 'rows' es el tamaño de la matriz de salida (definida previamente)
    D_voronoi_interp = zeros(rows);
    x_grid = linspace(-L/2, L/2, rows);
    y_grid = linspace(-L/2, L/2, rows);
    [YY, XX] = meshgrid(y_grid, x_grid);
    % Define el polígono rectangular de la zona de cobertura
    rect = polyshape([-L/2 L/2 L/2 -L/2], [-L/2 -L/2 L/2 L/2]);
    
    % --- 5. Rellenar la matriz de interpolación basándose en cada celda del Voronoi ---
    for i = 1:length(C)
        % Solo considerar celdas válidas (sin infinitos o NANs)
        if all(C{i} > 0)
            vx = V(C{i}, 1);  % Coordenadas x de los vértices de la celda
            vy = V(C{i}, 2);  % Coordenadas y de los vértices
            if any(isnan(vx)) || any(isnan(vy)) || any(isinf(vx)) || any(isinf(vy))
                continue;
            end
            % Crear el polígono Voronoi y recortarlo al rectángulo de la zona
            poly = polyshape(vx, vy);
            poly_clipped = intersect(poly, rect);
            if isempty(poly_clipped.Vertices)
                continue;
            end
            % Inicialmente se asume que la celda no detectó (valor 1)
            valor = 1;
            % Si el nodo detectó señal y NO fue marcado por el proceso de interpolación,
            % se fuerza el valor 0
            if decision_nodos(i) == 1 && ~marcar_blanco(i)
                valor = 0;
                  elseif nodos_aislados(i)
                valor = 1; % también zona gris
            end
            
            % Determinar las posiciones dentro del polígono recortado
            in_vec = isinterior(poly_clipped, XX(:), YY(:));
            in = reshape(in_vec, size(XX));
            % Asignar el valor a esa región de la malla
            D_voronoi_interp(in) = valor;
        end
    end
    
    % --- 6. Comparar la matriz interpolada (estimada) con la matriz real de SGH ---
    % Se asume que 'actual_SHG' es una matriz (rows x rows x runs) que indica los SGH reales
    xor_voronoi = not(xor(D_voronoi_interp, actual_SHG(:,:,run)));
    % Se calcula la tasa de acierto del SGH redondeado a dos decimales
    Frac_SHG_correct(run) = round(mean(xor_voronoi(:)), 2);
end


disp('******************************************************************');
disp( ['Parameter value: ', num2str(Parameter(par))]); % Display the parameter varied
disp( ['Measured SNR: ', num2str(10*log10(mean(mean(snr)))),'  dB']);
disp( ['Mean of SHGDR: ', num2str(mean(Frac_SHG_correct))]);
disp( ['Variance of SHGDR: ', num2str(var(Frac_SHG_correct))]);
disp( ['Min. local Pfa: ', num2str(min(Pfac)), '    Max. local Pfa: ', num2str(max(Pfac))]);
disp( ['Min. local Pd: ', num2str(min(Pdc)), '    Max. local Pd: ', num2str(max(Pdc))]);
disp('******************************************************************');

mean_SHGDR(par) = mean(Frac_SHG_correct);
var_SHGDR(par) = var(Frac_SHG_correct);
min_SHGDR(par) = min(Frac_SHG_correct);
max_SHGDR(par) = max(Frac_SHG_correct);

end % END OF PARAMETER (par) LOOP

%% Plot of the results
figure(1);
% Calculate the min and max differences
lowerError = min_SHGDR; % Difference to the min
upperError = max_SHGDR; % Difference to the max
errorbar(Parameter, mean_SHGDR, mean_SHGDR-lowerError, mean_SHGDR-upperError, 'bo-', ...
         'LineWidth', 1.5, 'MarkerFaceColor', 'w'); % Mean with min/max as error bars
hold on;
plot(Parameter, var_SHGDR, 'rs-', 'LineWidth', 1.5, 'MarkerFaceColor', 'w');
for i = 1:length(Parameter)
    text(Parameter(i), var_SHGDR(i)+0.02, sprintf('%.4f', var_SHGDR(i)), ... % Shift text downward
         'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
         'FontSize', 10, 'Color', 'r'); % Adjust appearance
end
grid on; axis auto;
ylabel('SHG detection rate');

if Flag == 1; xlabel('Total number of SSIoTs, {\it{N}}_1'); end
if Flag == 2; xlabel('M, {\it{M}}_1'); end
if Flag == 3; xlabel('SNR in dB'); end
if Flag == 4; xlabel('Numero de muestras'); end
if Flag == 5; xlabel('Factor em dtresh'); end

legend('Mean, Max and Min','Variance','Location','Best');

text(0.2, 0.3, {['$N_1=' num2str(N1) '$'], ...
               ['$m=' num2str(M) '$'], ...
                ['$s=' num2str(SNR) '$']}, ...
     'Interpreter', 'latex', ...
     'Units', 'normalized', ...
     'Color', 'blue', ...
     'FontSize', 12, ...
     'BackgroundColor', 'white');

hold off;


%% Exporting figure with tight borders and desired position and size
f = gcf; % Get the current figure handle
currentPos = f.Position; % Get the current position
f.Position = [currentPos(1), currentPos(2), 500, 300]; % Update only width and height
% folderPath = 'C:\Figures';
% fileName = 'SHGDR_versus_xxx.png';
% exportgraphics(gcf, fullfile(folderPath, fileName), 'Resolution', 600, 'ContentType', 'image');


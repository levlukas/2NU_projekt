%% ===== 2NU, semestralni projekt =====
% Zadani:
% Navrhnete a implementujte obdobu Gaussovy eliminacni metody, ktera
% resi soustavz s tridiagonalnimi maticemi.
%
% Autor:
% Lukas Lev, 2566660

%% Deklarace
n = 10;  % pocet rovnic
A = zeros(n);  % vstupni matice
b = zeros(n,1);


%% Funkce pro plneni diagonaly
% temp: naplneni diagonaly jen jednim cislem
function A = fill_diagonal(A,fill_low,fill_mid,fill_high,symmetric)
    % ARGUMENTY:
    %   A         ... vstupni matice
    %   fill_low  ... cislo / vektor pro vyplneni prvni vedlejsi diagonaly 
    %                 pod hlavni
    %   fill_mid  ... cislo / vektor pro vyplneni hlavni diagonaly
    %   fill_high ... cislo / vektor pro vyplneni prvni vedlejsi diagonaly 
    %                 nad hlavni
    %   symmetric ... nepovinne, pro naplneni symetricky
    
    % kontrola vstupu
    if ~ismatrix(A) && not(size(A,1) == size(1,A))
        error('Error: Input must be a square matrix.\n');
    end

    % pokud je vstupem cislo, napln nim vektor pro modularni naplneni
    % matice
    if ~isvector(fill_mid) & isnumeric(fill_mid)
        fill_mid = fill_mid(:);
    end
    if ~isvector(fill_low) & isnumeric(fill_low)
        fill_low = fill_low(:);
    end
    if ~isvector(fill_high) & isnumeric(fill_high)
        fill_high = fill_high(:); 
    end
    if not(isvector(fill_high) || isvector(fill_low) || isvector(fill_mid))
        error('Error: Input must be numeric.\n');
    end

    % kontrola symetrie
    if nargin < 5
        symmetric = false;
    end

    % naplneni matice
    A = zeros(size(A));
    for i = 1:size(A)
        if symmetric
            idx = min(i, size(A) - i + 1);  % zaruceni symetrie
        else
            idx = i;
        end
        A(i, i) = fill_mid(mod(idx-1, length(fill_mid)) + 1);
        if i > 1
            A(i, i-1) = fill_low(mod(idx-2, length(fill_low)) + 1);
        end
        if i < size(A)
            A(i, i+1) = fill_high(mod(idx-1, length(fill_high)) + 1);
        end
    end
end


%% Funkce pro naplneni vektoru prave strany
function b = fill_rhsvector(b,fill,symmetric)
    % ARGUMENTY:
    %   b         ... vstupni vektor
    %   fill      ... cislo / vektor pro vyplneni vektoru prave strany
    %   symmetric ... nepovinne, pro naplneni symetricky

    % kontrola vstupu
    if ~isvector(b)
        error('Error: Input must be a vector.\n');
    end
    if ~isvector(fill) & isnumeric(fill)
        fill = fill(:);
    elseif ~isvector(fill)
        error('Error: Input must be numeric.\n');
    end

    % kontrola symetrie
    if nargin < 3
        symmetric = false;
    end

    % naplneni vektoru
    for i = 1:length(b)
        if symmetric
                idx = min(i, length(b)-i+1);  % zaruceni symetri
        else
            idx = i;
        end
        b(i) = fill(mod(idx-1, length(fill)) + 1);
    end
end


% Funkce pro reseni systemu
function x = solve_tridiagonal(A, d)

    n = length(d);
    x = zeros(n,1);

    % rozpoznani diagonal
    diag_mid = zeros(n,1);
    diag_low = zeros(n-1,1);
    diag_high = zeros(n-1,1);
    for i = 1:n
        diag_mid(i) = A(i,i);
        if i > 1
            diag_low(i) = A(i,i-1);
        end
        if i < n
            diag_high(i) = A(i,i+1);
        end
    end

    % dopredna eliminace
    for i = 2:n
        if diag_mid(i-1) == 0  % zpetna kontrola pivotu
            error('System is not solvable: zero pivot at row %d', i-1);
        end
        m = diag_low(i-1) / diag_mid(i-1);
        diag_mid(i) = diag_mid(i) - m * diag_high(i-1);  % aktualizace pivotu
        d(i) = d(i) - m * d(i-1);  % aktualizace vektoru p. strany
    end

    % kontrola pivotu pro resitelnost
    if diag_mid(n) == 0
        error('System is not solvable: zero pivot at last row.');
    end

    % zpetna substituce
    x(n) = d(n) / diag_mid(n);
    for i = n-1:-1:1
        if diag_mid(i) == 0
            error('System is not solvable: zero pivot at row %d', i);
        end
        x(i) = (d(i) - diag_high(i) * x(i+1)) / diag_mid(i);
    end
end


%% main
A = fill_diagonal(A,-1,[3;4],-1);
b = fill_rhsvector(b,[2;1],true);

disp(A);

disp(solve_tridiagonal(A,b));

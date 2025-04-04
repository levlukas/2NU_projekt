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


%% Funkce pro plneni diagonaly
% temp: naplneni diagonaly jen jednim cislem
function A = fill_diagonal(A,fill_mid,fill_low,fill_high)
    % A ... vstupni matice
    % fill_mid ... cislo / vektor pro vyplneni hlavni diagonaly
    % fill_low ... cislo / vektor pro vyplneni prvni vedlejsi diagonaly pod hlavni
    % fill_high ... cislo / vektor pro vyplneni prvni vedlejsi diagonaly nad hlavni
    if ~isvector(fill_mid)
        % make into vector
    elseif ~isvector(fill_low)
    elseif ~isvector(fill_high)
    end
    for i = 0:size(A)
        A(i,i) = fill_mid(i);
        A(i+1,i) = fill_high(i);
        A(i-1,i) = fill_low(i);
    end
end


%% main
A = fill_diagonal(A,3);

disp(A);

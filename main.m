%% ===== 2NU, semestralni projekt =====
% Zadani:
% Navrhnete a implementujte obdobu Gaussovy eliminacni metody, ktera
% resi soustavz s tridiagonalnimi maticemi.
%
% Autor:
% Lukas Lev, 2566660

%% Deklarace
n = 100;  % pocet rovnic
A = zeros(n);  % vstupni matice


%% Funkce pro plneni diagonaly
function A = fill_diagonal(A,fill,diag_y)
    % A ... vstupni matice
    % fill ... cislo pro vyplneni diagonaly
    % diag_y ... souradnice diagonaly, jez ma byt vyplnena
    row = 0;
    while diag_y + row <= size(A)
        for row = 1:size(A)
            for col = 1:size(A)
                A(row,col) = 
            end
        end
    end
end

% --- LABORATION 1 ---
% @author Viola S�derlund
% @version 2020-03-21

% 1. Olinj�r skal�r ekvation

% Problem: Best�m samtliga r�tter till ekvationen:
% x^2 - 3sin(3x + 2) - 1 = 0

f = @(x) (x^2 - 3*sin(3*x + 2) - 1);
Df = @(x) (2*x - 9*cos(3*x + 2));

% a. Ritar graften till funktionen.

hold on
    range = [-0.8, 0.5];
    fplot(f, range)
    fplot(@(x) 0, range)
hold off

% b. Unders�ker vilka av funktionens r�tter som kan best�mmas av
%    fixpunktsiterationen next_x(x*) = x*, och ber�knar sedan r�tterna.

disp('------------------------------');
disp('--- b. Fixpunktiterationer ---');
disp('------------------------------');

next_x = @(x) (f(x) / 9 + x);
Dnext_x = @(x) (Df(x) / 9 + 1);

    % Fall 1
    disp('--- Fall 1 ---');
    
    xn = -0.6;
    x0 = xn
    
    [convergens, solution] = get_convergens(next_x, xn);
    
    if convergens
        if Dnext_x(solution) == 0
            disp('----> Iterationen konvergerar kvadratiskt.');
        else 
            disp('----> Iterationen konvergerar inte kvadratiskt.');
        end
    end
        
    % Fall 2  
    disp('--- Fall 2 ---');
    
    xn = 0.48;
    xo = xn
    
    [convergens, solution] = get_convergens(next_x, xn);
    
    if convergens
        if Dnext_x(solution) == 0
            disp('----> Iterationen konvergerar kvadratiskt.');
        else 
            disp('----> Iterationen konvergerar inte kvadratiskt.');
        end
    end
    
    disp('--> Correction: Andra nollv�rdet �r en falsk l�sning! Endast f�rsta nollv�rdet �r sann.')

% c. Ber�knar r�tterna med Newtons metod.

disp('------------------------');
disp('--- b. Newtons metod ---');
disp('------------------------');

next_x = @(x) (x - f(x) / Df(x));

    % Fall 1
    disp('--- Fall 1 ---');
    
    xn = -0.73;
    x0 = xn
    
    solution = calculate_solution(next_x, xn, next_x(xn));
    
    if Df(solution) ~= 0
        disp('----> Iterationen konvergerar kvadratiskt.');
    else 
        disp('----> Iterationen konvergerar inte kvadratiskt.');
    end
        
    % Fall 2
    disp('--- Fall 2 ---');
    
    xn = 0.47;
    x0 = xn
    
    solution = calculate_solution(next_x, xn, next_x(xn));
    
    if Df(solution) ~= 0
        disp('----> Iterationen konvergerar kvadratiskt.');
    else 
        disp('----> Iterationen konvergerar inte kvadratiskt.');
    end
    
    % Comparison
    
    disp('Fixpunktsiterationen (b) konvergerade l�ngsammare �n Newtons funktion (c), d� Newtons metod har en h�gre konvergensordning.');
    
% d. Beskriver Newtons metods konvergensordning.

disp('-------------------------------------------');
disp('--- d. Newtons metods konvergensordning ---');
disp('-------------------------------------------');

disp('Newtons metod har konvergensordning 2, d� fixpunktiterationen konvergerar linj�rt. Observera fall 1. Fixpunktsiterationen itererar 19 g�nger, medan Newton itererar 4 g�nger.');

fixed_point_iterations = 19
newton_method_iterations = 4

newton_method_convergence_order = floor(log(fixed_point_iterations) / log(newton_method_iterations))

function [convergens, solution] = get_convergens(next_x, xn)
       
    % Check convergens
    
    xn_incremented = next_x(xn);
        
    if abs(xn_incremented) > 1
        disp('Fixpunksiterationen divergerar runt xn = x0 ? x*.');
            
        convergens = false;
    elseif abs(xn_incremented) < 0.1
        disp('Fixpunktsiterationen konvergerar snabbt runt xn = x0 ? x*.');
            
        convergens = true;
    else
        disp('Fixpunktsiterationen konvergerar runt xn = x0 ? x*.')

        convergens = true;
    end
    
    % Calculate solution
    
    num_iterations = -1;
    solution = Inf;
        
    if convergens
        solution = calculate_solution(next_x, xn, xn_incremented);
    end
end

function solution = calculate_solution(next_x, xn, xn_incremented)
    
    % Calculate solution
    
    iterations = xn;
    solution = Inf;
        
    while xn_incremented ~= xn && abs(xn_incremented - xn) > eps
        xn = xn_incremented;
        xn_incremented = next_x(xn);
        
        iterations = [ iterations; xn ];
    end
    
    solution = xn;
    
    iterations
    num_iterations = length(iterations)
       
    solution
end
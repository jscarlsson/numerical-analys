% --- LABORATION 2.4bc ---
% @author Jakob Carlsson & whoever wrote the skeleton...
% @version 2020-04-19

function f=f_robotarm(z0,t,ts1,ts2,alpha,beta,gamma,omega)
    %% Computes the right-hand side of the robot differential equation
    %
    % * z0 is a vector containing [theta1(t);theta1prime(t);theta2(t);theta2prime(t)]
    % * t is the current time
    % ts1,ts2 is the two wanted angles.
    % 
    % The return value f is the right-hand side of the differential equation
    % (with four variables)
    % 

    % detta �r ist�llet f�r ts1 och ts2, men ok fair enough, vi k�r den och
    % sl�nger med dess output h�r ist�llet...
    theta_star_eh = get_theta(1.3,1.3);
    theta_star = theta_star_eh(:,1); % ta bara f�rsta kolumnen; detta borde inte vara n�dv�ndigt om man g�r r�tt men... eh.

    %initialisera f
    f=zeros(size(z0));
    % s�tt f, dvs systemet
    f(1)=z0(2);
    f(2)=-alpha*(z0(1)-theta_star(1) - gamma*z0(2) + beta*sin(omega*t));
    f(3)=z0(4);
    f(4)=-alpha*(z0(3)-theta_star(2) - gamma*(z0(4) + abs(z0(2))) + beta*sin(omega*t));
end

%forwards Euler, and also "animate" it by plotting it every now and then.
function animate_feuler(f)
    % r�kna ut tv� vinklar som en funktion av tiden ( = i i en loop)
    % vi har startv�rden, och vi har f fr�n f_robotarm som r�knar ut hf
    % (typ)



    %direkt fr�n slides:
    
end

function feuler(yprim, start, h)
    % given y' = yprim and y(0) = start, find y (a function)
    % we do this using forwards Euler with the step length h
    
    f  = yprim;
    y = start;
    for t=0:h:2-h
        y = y + h*f(t,y);
    end
end





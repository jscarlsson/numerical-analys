% --- LABORATION 2.4bc ---
% @author Jakob Carlsson & whoever wrote the skeleton...
% @version 2020-04-19

% this is kind of annoying, but to run you need to do:
% >> z0 = [pi/2;0;pi/6;0];
% >> f_robotarm(z0, 2, 0.5, 4, pi*3)

%I got rid of the parameters: ,t,ts1,ts2
function f=f_robotarm(z0,alpha,beta,gamma,omega)
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
    theta_star_eh = get_theta(1.3,1.3, false);
    theta_star = theta_star_eh(:,1); % ta bara f�rsta kolumnen; detta borde inte vara n�dv�ndigt om man g�r r�tt men... eh.

    %initialisera f
%     f=zeros(size(z0));
    % s�tt f, dvs systemet
    f{1}= @(t, z0) z0(2);
    f{2}= @(t, z0) -alpha*(z0(1)-theta_star(1)) - gamma*z0(2) + beta*sin(omega*t);
    f{3}= @(t, z0) z0(4);
    f{4}= @(t, z0) -alpha*(z0(3)-theta_star(2)) - gamma*(z0(4) + abs(z0(2))) + beta*sin(omega*t);
    
    h = 0.1;
    for t=0:h:15-h
        for i = 1:length(z0) %really we should be calling z0 just z at this point but whatever...
            z0(i) = z0(i) + h*f{i}(t,z0);
        end
        plot_robotarm([z0(1), z0(3)]); %OK SO, with 1 and 3, this gives the same result as time.m, but from the theory I THINK that it should actually be 2 and 4. But that looks worse so probably I'm wrong.
    end
end

% function feuler(yprim, start, h)
%     % given y' = yprim and y(0) = start, find y (a function)
%     % we do this using forwards Euler with the step length h
%     
%     f  = yprim;
%     y = start;
%     for t=0:h:2-h
%         y = y + h*f(t,y);
%     end
% end

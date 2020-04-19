% --- LABORATION 2.4bc ---
% @author Jakob Carlsson & whoever wrote the skeleton...
% @version 2020-04-18



% Funktionerna �r inte korrekta, man m�ste s�tta ihop dem p� korrekt
% s�tt, men jag tror detta �r r�tt byggstenar. Och man m�ste definiera u,
% som �r tv� olika (ish) i F_0 och F_1.


% eftersom Matlab b�rjar p� 1 m�ste vi ha u_1 och u_2
% f�rsta diffekvationen som system av f�rsta ordningens
F_0 = @(t, u) [u(2); -alpha*(u(1) - theta_star) - gamma*u(2) + beta*sin(omega*t)];

% andra diffek
F_1 = @(t, u) [u(2); -alpha*(u(1)-theta_star) - gamma*(u(2) + abs('F�RRA u(2)')) + beta*sin(omega*t)];

F = [F_0; F_1]; %ungef�r s� �r tanken, fast med parametrar d�...



% och just det.. theta_star �r vinklarna fr�n get_theta




% nu n�r jag kollar p� kommentaren f�r funktionen s� tror jag att jag vet
% hur man ska ta det jag skrev och faktiskt g�ra det anv�ndbart...
% Jag har lite fysiska anteckningar ocks� som jag inte vet hur bra de
% f�rmedlades h�r...


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


f=zeros(size(z0));
f(1)=z0(2);
f(2)=-alpha*(z0(1)-.... 
f(3)=z0(4);
f(4)=-alpha*(z0(3)-.... 



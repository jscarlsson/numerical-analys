% --- LABORATION 1 ---
% @author Viola S�derlund
% @version 2020-03-23

% 3. Olinj�rt ekvationsystem

A = struct('x', 93, 'y', 63, 'r', 55.1);
B = struct('x', 6, 'y', 16, 'r', 46.2);

% Newton's method

f = @(x, p) ((x(1) - p.x)^2 + (x(2) - p.y)^2 - p.r^2);
F = @(x) [ f(x, A) f(x, B) ];

f_x = @(x, p_x) (2*(x - p_x));
f_y = @(y, p_y) (2*(y - p_y));
G = @(x, p) [ f_x(x(1), p.x) f_y(x(2), p.y) ];
J = @(x) [ G(x, A); G(x, B) ];

xn = [ 34 54 ];
diff = [ -Inf, Inf ];

x0 = xn

while true
    diff = -F(xn)/J(xn) % J(xn)*sn = -F(xn), x(n+1) = xn + sn
    xn = xn + diff;
    
    if length(diff) <= 1
        break
    end
end

P_1 = struct('x', xn(1), 'y', xn(3))
P_2 = struct('x', xn(2), 'y', xn(4))

% Plotting the circles

hold on
    radians = 0:pi/50:2 * pi;
    x_unit = @(radius, offset_x) (radius*cos(radians) + offset_x);
    y_unit = @(radius, offset_y) (radius*sin(radians) + offset_y);
    
    radius = A.r; 
    plot(x_unit(radius, A.x), y_unit(radius, A.y));
    
    radius = B.r;
    plot(x_unit(radius, B.x), y_unit(radius, B.y));
    
    plot([ A.x B.x ], [ A.y B.y ], '.');
    plot([ P_1.x, P_2.x ], [ P_1.y, P_2.y ], '.');
hold off
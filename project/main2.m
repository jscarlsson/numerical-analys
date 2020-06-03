% --- PROJEKT KRETSEN ---
% @author Viola S�derlund
% @version 2020-05-11

load constants.mat;

% --- PLOT RESCALED CURRENT CURVE ---

h = 0.000001;

figure('Name', 'Omskalade strömkurvor', 'NumberTitle', 'off');

for i = 1:length(U0)
    [t, v] = vf(U0(i), F, h, 2);
    
    plot(t, v);
    hold on;
    plot(t, zeros(1, length(t)), '--', 'HandleVisibility', 'off'); % x-axis
    hold on;
end

legend('220 V', '1500 V', '2300 V', 'Location', 'northeast');

load constants.mat;

% --- ASSESS FOURIER SERIES KOEFFICIENTS ---

h = 0.000001;
num_k = 10;
num_p = 2;

a = [];

for i = 1:length(U0)
    [t, v] = vf(U0(i), F, h, num_p);
    [k, a_] = trf(t, v, num_k, h);
    
    a = [a a_];
    
    figure(U0(i));
    
    plot(t, v);
    hold on;
    plot(t, S(t, a_(1:3)))
    hold on;
    plot(t, S(t, a_));
    
    legend('Str�mkurvan (v)', 'Fourierutvecklingen (3)', 'Fourierutvecklingen (10)', 'Location', 'northeast');
end

disp('Series koefficients {a_k} for each U0:');
U0
a





load constants.mat;
% --- CONFIRM DECREASING WAVEFORMS ---

h = 0.000001;
num_k = 10;

odd = 1:2:num_k;

figure('Name', 'Semilogy-plotter över Fourierkoefficienterna', 'NumberTitle', 'off');

for i = 1:length(U0)
    [t, v] = vf(U0(i), F, h, 2);
    [k, a] = trf(t, v, num_k, h);
    
    semilogy(k(odd), abs(a(odd)));
    hold on;
end

legend('220 V', '1500 V', '2300 V', 'Location', 'northeast');





load constants.mat;
load mysterysound1.mat;

% --- INVESTIGATE MYSTERY SOUND ---

% Plot sound:

p = 2*pi;

h = p / length(v);
x = 0:h:p-h;
 

figure('Name', 'Ljudfilen', 'NumberTitle', 'off');


plot(x, v);
hold on;
plot(x, zeros(1, length(x)), '--', 'HandleVisibility', 'off'); % x-axis

legend('mysterysound1', 'Location', 'northeast');




load constants.mat;
load mysterysound1.mat;

% --- INVESTIGATE MYSTERY SOUND ---

h = 0.000001;
num_p = 2;


figure('Name', 'Mystery', 'NumberTitle', 'off');




% Calculate initial voltage U0* = U0_ of mystery sound:

% A graph of balanced current curves shows that
% the difference between the curves are the largest at x = pi/4.

% Results from 'fourier_mystery_sound_koefficients' shows that 
% 1500 V < U0* < 2300 V, as the koefficients decrement indicates
% the form of the current curve.

% Finding U0* with secant method based on the y-value at x = pi/4.

tgt = v(floor(length(v)/6)); % target y-value
eps = 10^-4;

xn_left = 1500; % = U0_left
xn_right = 2300; % = U0_right

last_approximation = Inf;
while true
    xn_new = next_xn(xn_left, xn_right, F, h, tgt, num_p);
    approximation = get_diff(xn_new, F, h, tgt, num_p);
   
    if approximation <= eps
        last_approximation = xn_new;
        break
    elseif approximation < 0
        xn_left = xn_new;
    else % approximation > 0
        xn_right = xn_new;
    end
        
    last_approximation = approximation;
end

U_ = last_approximation

% Confirm U*.

p = 2*pi;

h_ = p / length(v);
x = 0:h_:p-h_;

plot(x, v);

v_ = v;

hold on;

[t, v] = vf(U_, F, h, num_p);
plot(t, v);

hold on;
plot(x, zeros(1, length(x)), '--', 'HandleVisibility', 'off'); % x-axis

legend('Mystery sound', 'U0* V', 'Location', 'northeast');

% Due to the design of the method breakpoint, the error is at E(x = pi/4) <= eps = 10^-5
% However, the plot shows that the points do not form a smooth curve,
% therefore,

U_diff = [];
    
i = 1;
for xi = 1:length(x)
    target_t = x(xi);
        
    nearest_t = t(i);
    while t(i+1) < target_t
        i = i + 1;
        nearest_t = t(i);
    end
        
    if abs(nearest_t - target_t) > abs(t(i+1) - target_t)
        i = i + 1;
    end
        
    U_diff = [U_diff v(i)-v_(xi)];
end

% The largest error in this solution is,
U_diff_max = max(U_diff)

% Get difference between approximated value and target value.
function d = get_diff(U0, F, h, tgt, num_p)
    [t, v] = vf(U0, F, h, num_p);
    
    d = v(floor(length(v)/6)) - tgt;
end

function xn = next_xn(xn_left, xn_right, F, h, tgt, num_p)
    diff_left = get_diff(xn_left, F, h, tgt, num_p);
    diff_right = get_diff(xn_right, F, h, tgt, num_p);

    xn = xn_left - diff_left * (xn_left - xn_right) / (diff_left - diff_right);
end






% Fourier serie
function v = S(t, a)
    v = zeros(1, length(t));
    
    for i = 1:length(t)
        v(i) = 0;
        
        for k = 1:length(a)
            v(i) = v(i) + a(k)*sin(k*t(i));
        end
    end
end

% Trapezoidal rule
function I = T(h, k, t, v)
    f = @(i) v(i) * sin(k*t(i));

    I = f(1)/2;

    for i = 1:length(t)-1
        I = I + f(i);
    end
    
    I = h*(I + f(length(t))/2);
end
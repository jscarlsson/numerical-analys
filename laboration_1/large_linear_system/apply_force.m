load eiffel4.mat
style = "rand"; % "det" or "rand"
show = false; % whether to show a plot of the tower
numRuns = 20; % how many runs to take the average time of

n = length(xnod) %it's useful to show the size; double it to get N
forcepoint = n-50; %this is just to make it interesting, could be anything
b = choose_b(n, forcepoint, style); %forcepoint only relevant if "det"

disp('tid f�r inv(A)*b :');
for i=1:numRuns; tic; x = inv(A)*b; time = toc; end;
disp(time/numRuns);
disp('tid f�r A\b :');
for i=1:numRuns; tic; x = A\b; time = toc; end;
disp(time/numRuns);

if show
    xbel = xnod + x(1:2:end); ybel = ynod + x(2:2:end);
    hold on %clearly I don't get how hold works. This at least does what I want
        trussplot(xnod, ynod, bars);
        hold on
        trussplot(xbel, ybel, bars);
        if style=="det"
            hold on
            plot(xnod(forcepoint), ynod(forcepoint), '*');
            %hold on
            plot(xbel(forcepoint), ybel(forcepoint), '*');
        end
    hold off
end

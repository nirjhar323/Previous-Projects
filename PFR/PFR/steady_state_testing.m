% Steady state testing

load steady_state_testing.txt;

[m,n] = size(steady_state_testing);

x = linspace(1,m,m);

for i = 1:n
    
    plot(x,steady_state_testing(:,i));
    hold on;
end
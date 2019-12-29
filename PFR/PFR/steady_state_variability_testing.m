%Steady State Testing Comparison

load steady_state_testing.txt;

[m,n] = size(steady_state_testing);

load steady_state_testing2.txt;

[p,q] = size(steady_state_testing2);

x = linspace(1,m,m);
z = linspace(1,p,p);

 
    figure;
    plot(x,steady_state_testing(:,1));
    hold on;
    plot(z,steady_state_testing2(:,1));
    
    figure;
    plot(x,steady_state_testing(:,2));
    hold on;
    plot(z,steady_state_testing2(:,2));
    
    
    figure;
    plot(x,steady_state_testing(:,3));
    hold on;
    plot(z,steady_state_testing2(:,3));
    
    
    figure;
    plot(x,steady_state_testing(:,4));
    hold on;
    plot(z,steady_state_testing2(:,4));
    
    
    figure;
    plot(x,steady_state_testing(:,5));
    hold on;
    plot(z,steady_state_testing2(:,5));
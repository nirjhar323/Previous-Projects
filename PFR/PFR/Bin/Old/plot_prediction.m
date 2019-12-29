load data.txt;
for i = 1:501
    x(i) = 0.01*i - 0.01;
end

    [m,n] = size(data);
    
    
    for a = 1:n
        for b = 1:number
    plot(x,data(:,a));
    hold on;
    end
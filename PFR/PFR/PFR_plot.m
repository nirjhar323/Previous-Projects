%% Plot PFR data
%% Define plot settings here - Can support upto 100 CSTR plugs with 100 sets of parameters(i.e no_of_plots) in the text files
no_of_plugs = 30; %30;
no_of_plots = 80; %80;
line_width = 1.5;
edge_alpha = 0.01;


%%
load CH4.txt;
CH4 = CH4(:,1:no_of_plugs);

x = linspace(1,no_of_plugs,no_of_plugs);

figure; 
for a = 1:no_of_plots
patchline(x,CH4(a,:),'edgecolor','red','linewidth',line_width,'edgealpha',edge_alpha);
end
xlabel('No of Plugs');ylabel('Partial Pressure(bar)');title('CH4');
xlim([1 no_of_plugs]);
savefig('CH4.fig');
%%
load H2O.txt;

H2O = H2O(:,1:no_of_plugs);

x = linspace(1,no_of_plugs,no_of_plugs);

figure;
for a = 1:no_of_plots
patchline(x,H2O(a,:),'edgecolor','red','linewidth',line_width,'edgealpha',edge_alpha);
end
xlabel('No of Plugs');ylabel('Partial Pressure(bar)');title('H2O');
xlim([1 no_of_plugs]);
savefig('H2O.fig');
%%
load CO2.txt;


CO2 = CO2(:,1:no_of_plugs);

x = linspace(1,no_of_plugs,no_of_plugs);

figure;
for a = 1:no_of_plots
patchline(x,CO2(a,:),'edgecolor','red','linewidth',line_width,'edgealpha',edge_alpha);
end
xlabel('No of Plugs');ylabel('Partial Pressure(bar)');title('CO2');
xlim([1 no_of_plugs]);
savefig('CO2.fig');
%%

load H2.txt;

H2 = H2(:,1:no_of_plugs);

x = linspace(1,no_of_plugs,no_of_plugs);

figure;
for a = 1:no_of_plots
patchline(x,H2(a,:),'edgecolor','red','linewidth',line_width,'edgealpha',edge_alpha);
end
xlabel('No of Plugs');ylabel('Partial Pressure(bar)');title('H2');
xlim([1 no_of_plugs]);
savefig('H2.fig');
%%

load CO.txt;

CO = CO(:,1:no_of_plugs);

x = linspace(1,no_of_plugs,no_of_plugs);

figure
for a = 1:no_of_plots
patchline(x,CO(a,:),'edgecolor','red','linewidth',line_width,'edgealpha',edge_alpha);
end
xlabel('No of Plugs');ylabel('Partial Pressure(bar)');title('CO');
xlim([1 no_of_plugs]);
savefig('CO.fig');
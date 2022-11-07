% plotClCd.m

clearvars; close all; clc;

nrels814 = importdata('NRELs814.txt');
sg6040 = importdata('sg6040.txt');
S814_15_100Hz = importdata('S814_15_100Hz.txt');
naca63815 = importdata('naca63815.txt');

hfig = figure('Position',[100 100 800 400]);
hax = axes('Parent',hfig);
plot(nrels814(:,1),nrels814(:,2),'-r','DisplayName','S814 Lift')
hold on
plot(sg6040(:,1),sg6040(:,2),'-b','DisplayName','SG6040 Lift')
%plot(S814_15_100Hz(:,1),S814_15_100Hz(:,2),'-c','DisplayName','S814 v2 Lift')
%plot(naca63815(:,1),naca63815(:,2),'-k','DisplayName','NACA63815 Lift')
plot(nrels814(:,1),nrels814(:,3),'--r','DisplayName','S814 Drag')
plot(sg6040(:,1),sg6040(:,3),'--b','DisplayName','SG6040 Drag')
%plot(S814_15_100Hz(:,1),S814_15_100Hz(:,3),'--c','DisplayName','S814 v2 Drag')
%plot(naca63815(:,1),naca63815(:,3),'--k','DisplayName','NACA63815 Drag')
hold off
xlabel('AoA (deg)');
ylabel('Coefficient Value');
legend('Location','best')
export_fig(hfig,[pwd '\..\..\..\figures\liftDrag.png'],'-transparent','-m3');

hfig = figure('Position',[100 100 800 400]);
hax = axes('Parent',hfig);
plot(nrels814(:,1),nrels814(:,2)./nrels814(:,3),':r','Linewidth',2.0,'DisplayName','S814')
hold on
plot(sg6040(:,1),sg6040(:,2)./sg6040(:,3),':b','Linewidth',2.0,'DisplayName','SG6040')
%plot(S814_15_100Hz(:,1),S814_15_100Hz(:,2)./S814_15_100Hz(:,3),':c','Linewidth',2.0,'DisplayName','S814 v2')
%plot(naca63815(:,1),naca63815(:,2)./naca63815(:,3),':k','Linewidth',2.0,'DisplayName','NACA63815')
hold off
xlabel('AoA (deg)');
ylabel('Lift to Drag Ratio');
legend('Location','best')
export_fig(hfig,[pwd '\..\..\..\figures\liftToDrag.png'],'-transparent','-m3');
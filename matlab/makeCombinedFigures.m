% makeCombinedFigures.m

clearvars; close all; clc;

basepath = 'output/';
%baseGroup = 'MB_F1d5_S814'; endIndex = 45; Ylim = 3500;
baseGroup = 'MB_F1d5_SG6040_TSR2';  endIndex = 91; Ylim = 12000;
%baseGroup = 'MB_F1d5_SG6040'; endIndex = 71; Ylim = 8000;
%baseGroup = 'EastRiver_0_35'; 

% open the csv files and get data and maximums
pitch = [0 10];
for i=1:1:numel(pitch)
    simGroupDirName = [baseGroup '_p' char(string(pitch(i)))];
    T = readtable([basepath simGroupDirName '/' simGroupDirName '.xlsx']);
    netPwr(i,:) = T.netPower;
    towSpeed(i,:) = T.towspeed;
    noTowPwr(i) = T.netPower(1);
    % find ros that has max power
    [mxPwr(i),row] = max(T.netPower);
    towSpeedMaxPwr(i) = T.towspeed(row);
    inductionMaxPwr(i) = T.induction(row);
    thrustMaxPower(i) = T.thrust(row);
    percentOverNoTow(i) = T.powerGain(row)/noTowPwr(i)*100;
    percentOverNoTowPitch0(i) = (mxPwr(i) - noTowPwr(i))/noTowPwr(i)*100;
    percentOfBetzLimit(i) = mxPwr(i)/(16/27*0.5*998*1.5^3*pi)*100;
end
color = 'none';
hfig = figure
hfig.Color = color;
plot(-towSpeed(1,1:endIndex),netPwr(1,1:endIndex),'dk','LineWidth',2.0,'DisplayName','No Pitch');
hold on
plot(-towSpeed(2,1:endIndex),netPwr(2,1:endIndex),'ok','LineWidth',2.0,'DisplayName','10 Degrees');
%plot(-towSpeedMaxPwr(1),mxPwr(1),'dr','LineWidth',2.0,'DisplayName','');
%plot(-towSpeedMaxPwr(2),mxPwr(2),'or','LineWidth',2.0,'DisplayName','');
grid on
ylabel('Net Power (W)');
xlabel('Tow Speed (m/s)');
ax = gca;
ax.YLim = [-1000 Ylim];
ax.Color = color;
legend('Location','best');
export_fig(hfig,['figures\' baseGroup '_combinedFigureNetPower.png'],'-transparent','-m3');
saveas(hfig,['figures\' baseGroup '_combinedFigureNetPower.fig']);
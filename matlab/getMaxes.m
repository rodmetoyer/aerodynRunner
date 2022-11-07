% getMaxes.m

clearvars; close all; clc;

basepath = 'output/';
%baseGroup = 'MB_F1d5_S814';
%baseGroup = 'MB_F1d5_SG6040_TSR2';
baseGroup = 'MB_F1d5_SG6040';
%baseGroup = 'EastRiver_0_35';

% Get pitch 0 baseline data
%p0 = load([basepath baseGroup '_p0_towSpeed0/EastRiver_0_35_p0_towSpeed0.mat']);
%p0 = load([basepath baseGroup '_p0_towSpeed0/MB_F1d5_S814_p0_towSpeed0.mat']);
%p0 = load([basepath baseGroup '_p0_towSpeed0/MB_F1d5_SG6040_TSR2_p0_towSpeed0.mat']);
p0 = load([basepath baseGroup '_p0_towSpeed0/MB_F1d5_SG6040_p0_towSpeed0.mat']);

[baseMaxPower,indPwr] = max(p0.parsedResults.RtAeroPwr);

[hax, ~] = tight_subplot(1, 2, [.1 .1], [.1 .1], [.1 .1]);
hfig = hax(1).Parent;
hfig.Position = [100 100 900 400];
plot(hax(1),p0.parsedResults.RtTSR,p0.parsedResults.RtAeroPwr,'*k','LineWidth',2.0);
hold(hax(1),'on');
plot(hax(1),p0.parsedResults.RtTSR(indPwr),p0.parsedResults.RtAeroPwr(indPwr),'sr','MarkerSize',10.0,'LineWidth',2.0);
hold(hax(1),'off');
xlabel(hax(1),'TSR');
ylabel(hax(1),'Power (W)');
hax(1).Color = 'none';
plot(hax(2),p0.parsedResults.RtTSR,p0.parsedResults.RtAeroFxh,'*k','LineWidth',2.0);
hold(hax(2),'on');
plot(hax(2),p0.parsedResults.RtTSR(indPwr),p0.parsedResults.RtAeroFxh(indPwr),'sr','MarkerSize',10.0,'LineWidth',2.0);
hold(hax(2),'off');
xlabel(hax(2),'TSR');
ylabel(hax(2),'Thrust (N)');
grid(hax(1),'on'); grid(hax(2),'on');
hax(2).Color = 'none';
export_fig(hfig,['figures\' baseGroup 'exampleResults.png'],'-transparent','-m3');

% figure
% plot(p0.parsedResults.RtTSR,p0.parsedResults.RtAeroCp,'.k');

% open the csv files and get the maximums
pitch = 0:1:20;
for i=1:1:numel(pitch)
    simGroupDirName = [baseGroup '_p' char(string(pitch(i)))];
    T = readtable([basepath simGroupDirName '/' simGroupDirName '.xlsx']);
    noTowPwr(i) = T.netPower(1);
    % find ros that has max power
    [mxPwr(i),row] = max(T.netPower);
    towSpeedMaxPwr(i) = T.towspeed(row);
    inductionMaxPwr(i) = T.induction(row);
    thrustMaxPower(i) = T.thrust(row);
    percentOverNoTow(i) = T.powerGain(row)/noTowPwr(i)*100;
    percentOverNoTowPitch0(i) = (mxPwr(i) - baseMaxPower)/baseMaxPower*100;
    percentOfBetzLimit(i) = mxPwr(i)/(16/27*0.5*998*1.5^3*pi)*100;
end

hfig = figure
hfig.Color = 'none';
plot(pitch,percentOfBetzLimit,'*k','LineWidth',2.0);
xlabel('Pitch (deg)');
ylabel('Percent Betz Limit');
ax = gca;
ax.Color = 'none';
export_fig(hfig,['figures\' baseGroup 'percentOfBetzLimit.png'],'-transparent','-m3');

hfig = figure
hfig.Color = 'none';
plot(pitch,mxPwr,'*k','LineWidth',2.0);
xlabel('Pitch (deg)');
ylabel('Max Power (W)');
ax = gca;
ax.Color = 'none';
export_fig(hfig,['figures\' baseGroup 'maxPwr.png'],'-transparent','-m3');

hfig = figure
hfig.Color = 'none';
plot(pitch,noTowPwr,'*k','LineWidth',2.0);
xlabel('Pitch (deg)');
ylabel('Power When Stationary');
ax = gca;
ax.Color = 'none';
export_fig(hfig,['figures\' baseGroup 'noTowPwr.png'],'-transparent','-m3');

hfig = figure
hfig.Color = 'none';
plot(pitch,percentOverNoTow,'*k','LineWidth',2.0);
xlabel('Pitch (deg)');
ylabel('Percent Increase Over Stationary (%)');
ax = gca;
ax.Color = 'none';
export_fig(hfig,['figures\' baseGroup 'percentOverNoTow.png'],'-transparent','-m3');

hfig = figure
ax.Color = 'none';
plot(pitch,percentOverNoTowPitch0,'*k','LineWidth',2.0);
xlabel('Pitch (deg)');
ylabel('Percent Increase Over Stationary and No Pitch (%)');
ax = gca;
ax.Color = 'none';
export_fig(hfig,['figures\' baseGroup 'percentOverNoTowPitch0.png'],'-transparent','-m3');

hfig = figure
hfig.Color = 'none';
plot(pitch,inductionMaxPwr,'*k','LineWidth',2.0);
xlabel('Pitch (deg)');
ylabel('Induction Factor at Max Power');
ax = gca;
ax.Color = 'none';
export_fig(hfig,['figures\' baseGroup 'inductionMaxPwr.png'],'-transparent','-m3');

return

hfig = figure
hfig.Color = 'none';
plot(pitch,towSpeedMaxPwr,'.k','LineWidth',2.0,'MarkerSize',8.0);
xlabel('Pitch (deg)');
ylabel('Tow Speed at Max Power (m/s)');
ax = gca;
ax.Color = 'none';
export_fig(hfig,['figures\' baseGroup 'towSpeedMaxPwr.png'],'-transparent','-m3');









hfig = figure
hfig.Color = 'none';
plot(pitch,thrustMaxPower,'.k','LineWidth',2.0,'MarkerSize',8.0);
xlabel('Pitch (deg)');
ylabel('Thrust at Max Power (N)');
ax = gca;
ax.Color = 'none';
export_fig(hfig,['figures\' baseGroup 'thrustMaxPower.png'],'-transparent','-m3');

hfig = figure
hfig.Color = 'none';
plot(towSpeedMaxPwr,mxPwr,'.k','LineWidth',2.0,'MarkerSize',8.0);
xlabel('Tow Speed at Max Power (m/s)');
ylabel('Max Power (W)');
ax = gca;
ax.Color = 'none';
export_fig(hfig,['figures\' baseGroup 'towSpeedMaxPwrmxPwr.png'],'-transparent','-m3');
% computeTwist.m
% Computes the optimal twist for a blade.
% Option to call makeAirfoilDataFiles.m to make the data files necessary
% for making the blade in SolidWorks

% Rodney Metoyer 2018
% North Carolina State University
% This code is free to use and distribute

clearvars; close all;
savePlots = true;

%% Design parameters
numBlades = 3;
bladeRadius_m = 0.3556*0.5;
aspectRatio = 10;
bladeChord_m = bladeRadius_m/aspectRatio;%0.625*0.0254; % Inches times meters per inch
bladeDZfrac = 0.0; % Blade dead zone fraction (Ro/R - see Spera 2009)
numSections = 50; % Number of sections (whole number)
rtos = 2.0;    % This is the ratio of rotor radius to disturbed fluid
% stream length for computing optimal TSR. A good rule-of-thumb value is
% 2.0 (Ragheb and Ragheb 2011).
AoAopt_deg = 8.0; % Optimal angle of attack for the airfoil
pitch = 0;

%% Compute optimal TSR
% Using method described in (Ragheb and Ragheb 2011)
TSRarray = 0.1; %2*pi/numBlades*rtos;
name = ['S814_windTunnelRotor_g'];
%TSR = 1.0;
%TSRarray = [0.1 0.75 2.0 2*pi/numBlades*rtos];

%% Compute Twist
hfig = figure;
ax = axes('Parent',hfig);
hold(ax,'on');
markers = '*+dx';
for TSRitr = 1:1:numel(TSRarray)
    TSR = TSRarray(TSRitr);
    % Using method described in Ch. 5 of (Gasch and Twele 2012)
    locs_m = linspace(bladeRadius_m*bladeDZfrac,bladeRadius_m,numSections);
    twist_deg = atand(2/3*bladeRadius_m/TSR*1./locs_m)-AoAopt_deg + pitch;

    %% Make files
    if true % Change to false if you don't want to make files
        basefile = [pwd '\simulationInputFiles\S814_base.txt'];
        
        % convert to mm
        % bladeChord_mm = bladeChord_m*1000.0;
        %locs_mm = locs_m*1000.0;
        [x,y,z] = getCoords(basefile,name,bladeChord_m,locs_m,twist_deg);
    end

    %% Plots
       figure
       plot(ax,locs_m,twist_deg,markers(TSRitr),'DisplayName',['TSR = ' num2str(TSR,'%3.2f')]);
       xlabel(ax,'Location (m)'); ylabel(ax,'Twist Angle (deg)');
       hfig2 = figure;
       ax2 = axes('Parent',hfig2);
       plot3(ax2,x(:,1),y(:,1),z(:,1));
       hold on
       axis equal;
       for i=2:1:length(twist_deg)
           plot3(ax2,x(:,i),y(:,i),z(:,i));
       end
       view(-21,61);
       hold off
       ax2.Color = 'none';
       ax2.XTickLabel = '';
       ax2.YTickLabel = '';
       ax2.ZTickLabel = '';
    if savePlots
        if ~exist([pwd '\figures\'], 'dir')
            mkdir([pwd '\figures\']);
        end
       export_fig(hfig2,[pwd '\figures\TSR' num2str(TSR,'%3.1f') 'bladePlot.png'],'-transparent','-m3');
    end
    T = table(locs_m.',twist_deg.');
    writetable(T,[pwd '\figures\' name '.xlsx']);
%     locs_m,twist_deg
%     writematrix()
end
ax.Color = 'none';
ax.XLim = [0 round(bladeRadius_m*1.1,2)];
ax.YLim = [0 90];
legend(ax,'Color','w','Location','SouthWest');
if savePlots
    export_fig(hfig,[pwd '\figures\TSRtwistDistAll.png'],'-transparent','-m3');
end

function [x,y,z] = getCoords(basefile,name,c,locs,twistAngle_deg)
    % Makes data files for creating blades/wings in solidworks from airfoils at
    % spanwise stations

    % INPUT - The basefile for the airfoil that contains the x/c and y/c coords

    chord = c;
    centerLoc = 0.25;
    numSpanwiseLocs = length(twistAngle_deg);
    if length(locs) ~= numSpanwiseLocs
        error('Number of locations and twist angles must match');
    end

    A = importdata(basefile);

    %x = A(:,1)*chord;
    % Shifting back so that 1/4 chord is located at orgin
    A(:,1) = A(:,1)-centerLoc;
    x = A(:,1)*chord*cosd(twistAngle_deg)-A(:,2)*chord*sind(twistAngle_deg);
    %y = A(:,2)*chord;
    y = A(:,1)*chord*sind(twistAngle_deg)+A(:,2)*chord*cosd(twistAngle_deg);

%     pfold = [pwd '\airfoilData\' name];
%     if ~exist(pfold, 'dir')
%         mkdir(pfold);
%     end

    for i=1:1:numSpanwiseLocs
        z(:,i) = ones(length(x),1)*locs(i);
        B = [x(:,i) y(:,i) z(:,i)];
        if ~exist([pwd '\figures\'], 'dir')
            mkdir([pwd '\figures\']);
        end
        filename = [pwd '\figures\' name '_' num2str(i) '.txt']; 
        fid = fopen(filename,'w');
        fprintf(fid,'%12.10f %12.10f %12.10f\n',B');
        fclose(fid);
    end
end
% findBestStationaryMeterRadius.m
% Runs aerodyn analyses for the towed turbine project.
% Written by Rodney Metoyer
% All rights reserved


%% Clear the workspace
clearvars; close all; clc;

%% Set up simulations
flowspeed = 1.0; % This is the speed of the flow, not the relative flow speed seen by the turbine.
makeAllThePlots = false;

simGroupName = 'Looking for best still rotor'; simGroupDirName = 'bestStill';
sims(1).name = "still_pm5";
sims(2).name = "still_pm4"; 
sims(3).name = "still_pm3"; 
sims(4).name = "still_pm2"; 
sims(5).name = "still_pm1"; 
sims(6).name = "still_p0"; 
sims(7).name = "still_p1"; 
sims(8).name = "still_p2";
sims(9).name = "still_p3"; 
sims(10).name = "still_p4";
sims(11).name = "still_p5"; 

runAerodynSimulations = true;
parseAerodynResults = true;

if ~strcmp(pwd,"D:\School\Dissertation\TranslatingTurbine\newAerodyn\matlab")
    error('Looks like a previous run was interupted');
end

%% Run AeroDyn
if(runAerodynSimulations)
    currDir = cd([pwd '\simulationInputFiles']);
    for i=1:1:numel(sims)
        driverFileNames = [char(sims(i).name) '_driver.dvr'];
        command=['..\..\bin\AeroDyn_Driver_Win32 ' driverFileNames];
        status = dos(command);
        % move results into their own folder
        if ~exist(['..\output\' char(sims(i).name)],'dir')
            mkdir(['..\output\' char(sims(i).name)]);
        end
        disp(['Moving results for ' char(sims(i).name) ', please stand by']);
        [status,message,messageId] = movefile('*.out',['..\output\' char(sims(i).name)]);
        pause(1);
        if ~status
            disp(message);
        end
    end
    cd(currDir);
end

%% Process results
pathToOutputFolder = [pwd '\output'];
if(parseAerodynResults)    
    for i=1:1:numel(sims)
        disp(['Parsing results for ' char(sims(i).name) ', please stand by']);
        pathToFolder = [pathToOutputFolder '\' char(sims(i).name)];
        sims(i).results = getResults(pathToFolder);        
        for j=1:1:numel(sims(i).results)
            towspeed(j) = sims(i).results(j).fluidSpeed - flowspeed;
            pitch(j) = sims(i).results(j).pitch;
            for k=1:1:numel(sims(i).results(j).headers)
                if strcmp('RtTSR',sims(i).results(j).headers(k))
                    sims(i).RtTSR(j) = mean(sims(i).results(j).data(:,k));
                elseif strcmp('RtAeroPwr',sims(i).results(j).headers(k))
                    sims(i).RtAeroPwr(j) = mean(sims(i).results(j).data(:,k));
                elseif strcmp('RtAeroFxh',sims(i).results(j).headers(k))
                    sims(i).RtAeroFxh(j) = mean(sims(i).results(j).data(:,k));
                elseif strcmp('RtSpeed',sims(i).results(j).headers(k))
                    sims(i).RtSpeed(j) = mean(sims(i).results(j).data(:,k));
                elseif strcmp('RtAeroCp',sims(i).results(j).headers(k))
                    sims(i).RtAeroCp(j) = mean(sims(i).results(j).data(:,k));
                elseif strcmp('RtAeroCt',sims(i).results(j).headers(k))
                    sims(i).RtAeroCt(j) = mean(sims(i).results(j).data(:,k));
                end
                %sims(i).(sims(i).results(j).headers(k)) = sims(i).results(j).data(:,k);
            end
        end
        towspeed = unique(towspeed);
        if numel(towspeed)>1
            error('Should all be the same speed for every simulation case.');
        end
        pitch = unique(pitch);
        if numel(pitch)>1
            error('Should all be the same speed for every simulation case.');
        end
        sims(i).towspeed = towspeed;
        sims(i).pitch = pitch;
        parsedResults = sims(i);
        save([pathToFolder '\' char(sims(i).name) '.mat'],'parsedResults');
    end    
end

% If not parsing, assume we are loading
if ~parseAerodynResults
    for i=1:1:numel(sims)
        pathToFolder = [pathToOutputFolder '\' char(sims(i).name)];
        tempsims(i) = load([pathToFolder '\' char(sims(i).name) '.mat']);
    end
    clear sims;
    for i=1:1:numel(tempsims)
        sims(i) = tempsims(i).parsedResults;
    end    
    clear tempsims;
end

%% Compute max power, thrust at max power, and induction factor at max power
syms a
for i=1:1:numel(sims)
    for j=1:1:numel(sims(i).RtAeroCp)
        cp = sims(i).RtAeroCp(j);
        S = vpasolve(4*a*(1-a)^2 == cp, a);
        S = double(S);
        S = S(S<=1.0);
        if cp > 0
            S = S(S>=0.0);
        end
        if numel(S) > 1
            S = S(S<0.5);
        end
        if numel(S) ~= 1
            disp('empty');
        end
        sims(i).anum(j) = S;
    end
    % max power
    [sims(i).maxCp,ind] = max(sims(i).RtAeroCp);
    [sims(i).maxPower,ind2] = max(sims(i).RtAeroPwr);
    if ind ~= ind2
        error('This is a problem. They should be in the same place');
    end
    sims(i).thrustMaxCp = sims(i).RtAeroFxh(ind);
    sims(i).indFactMaxCp = sims(i).anum(ind);
    
    casename{i} = char(sims(i).name);
    power(i) = sims(i).maxPower;
    thrust(i) = sims(i).thrustMaxCp;
    induction(i) = sims(i).indFactMaxCp;
    towPower(i) = sims(i).thrustMaxCp*sims(i).towspeed;
    netPower(i) = sims(i).maxPower - sims(i).thrustMaxCp*sims(i).towspeed;
    % warning, the next line assumes still is the first simulation
    powerGain(i) = netPower(i) - sims(1).maxPower;
    towspeed(i) = sims(i).towspeed;
    pitch(i) = sims(i).pitch;
    
    disp([char(sims(i).name) ' Power: ' num2str(power(i),'%6.2f') ' (W)']);
    disp([char(sims(i).name) ' Thrust: ' num2str(thrust(i),'%6.2f') ' (N)']);
    disp([char(sims(i).name) ' Induction: ' num2str(induction(i),'%6.4f') ' (unitless)']);
    disp([char(sims(i).name) ' Tow Power: ' num2str(towPower(i),'%6.4f') ' (W)']);
    disp([char(sims(i).name) ' Net Power: ' num2str(netPower(i),'%6.4f') ' (W)']);
    disp([char(sims(i).name) ' Power Gain: ' num2str(powerGain(i),'%6.4f') ' (W)']);
    disp([char(sims(i).name) ' Power Gain %: ' num2str(powerGain(i)/sims(1).maxPower*100,'%6.4f') ' (W)']);
    disp(' ');
end
power = power.'; thrust = thrust.'; induction = induction.'; towPower = towPower.';
netPower = netPower.'; powerGain = powerGain.'; towspeed = towspeed.';
T = table(towspeed,power,thrust,induction,towPower,netPower,powerGain,'RowNames',casename);
% write this to file
writeDir = [pwd '\output\' simGroupDirName];
if ~exist(writeDir,'dir')
    mkdir(writeDir);
end
writetable(T,[writeDir '\' simGroupDirName '.xlsx'],'WriteRowNames',true);


%% Make plots
if ~exist([writeDir '\figures'],'dir')
    mkdir([writeDir '\figures']);
end
if makeAllThePlots
    for i=1:1:numel(sims)
        hfig = figure;
        ax = axes('Parent',hfig);
        plot(ax,sims(i).RtTSR,sims(i).RtAeroPwr,'*');
        xlabel(ax,'TSR'); ylabel(ax,'Power (W)');
        title(ax,sims(i).name,'Interpreter','none');
        ax.Color = 'none';
        export_fig(hfig,[writeDir '\figures\pwrVtsr.png'],'-transparent','-m3');

        hfig = figure;
        ax = axes('Parent',hfig);
        plot(ax,sims(i).RtTSR,sims(i).RtAeroCp,'*');
        xlabel(ax,'TSR'); ylabel(ax,'Power Coefficient (c_P)');
        title(ax,sims(i).name,'Interpreter','none');
        ax.Color = 'none';
        export_fig(hfig,[writeDir '\figures\CpVtsr.png'],'-transparent','-m3');

        hfig = figure;
        ax = axes('Parent',hfig);
        plot(ax,sims(i).RtTSR,sims(i).RtAeroFxh,'*');
        xlabel(ax,'TSR'); ylabel(ax,'Thrust (N)');
        title(ax,sims(i).name,'Interpreter','none');
        ax.Color = 'none';
        export_fig(hfig,[writeDir '\figures\ThrustVtsr.png'],'-transparent','-m3');

        hfig = figure;
        ax = axes('Parent',hfig);
        plot(ax,sims(i).RtTSR,sims(i).anum,'*');
        xlabel(ax,'TSR'); ylabel(ax,'Induction Factor');
        title(ax,sims(i).name,'Interpreter','none');    
        ax.Color = 'none';
        export_fig(hfig,[writeDir '\figures\InductionVtsr.png'],'-transparent','-m3');
    end
end

hfig = figure;
ax = axes('Parent',hfig);
plot(ax,pitch,netPower,'o');
grid(ax,'on');
xlabel(ax,'Pitch (deg)'); ylabel(ax,'Power (W)');
title(ax,'Power vs. Pitch for Stationary Turbine');
ax.Color = 'none';
export_fig(hfig,[writeDir '\figures\pitchVnetPower.png'],'-transparent','-m3');



% Get results from the results files and puts them into a struct
function results = getResults(pathToFolder)
    currDir = cd(pathToFolder);
    resultsList = dir('*.out');
    disp([pathToFolder ' has ' num2str(numel(resultsList),'%3i') ' output files.']);
    for i=1:1:numel(resultsList)
        outfile = [resultsList(i).name];
        % Every simualtion case has an output file.
        fid = fopen(outfile);      
        if fid == -1
            error(['Unable to open file: ' outfile]);
        end
        % Read file a line at a time
        lineOfText = fgetl(fid);
        lineNumber = 1;
        dataRow = 1;
        while ~isnumeric(lineOfText)            
            % TODO get case information
            if lineNumber == 5
                % Line 5 has case data
                caseData = strsplit(lineOfText,[" "],'CollapseDelimiters',true);
                results(i).caseNumber = str2double(caseData{2}(1:end-1));
                [~,spd] = strtok(caseData{3},'=');
                results(i).fluidSpeed = str2double(spd(2:end));
                clear spd;
                [~,spd] = strtok(caseData{7},'=');
                results(i).rotorSpeed = str2double(spd(2:end));
                clear spd;
                [~,spd] = strtok(caseData{9},'=');
                results(i).pitch = str2double(spd(2:end));
            elseif lineNumber == 7
                % Line 7 has the headers (TODO make this safer - check for a header string)
                headers = strsplit(lineOfText,[" ","\t"],'CollapseDelimiters',true);
                results(i).headers = string(headers(2:end-1));
            elseif lineNumber == 8
            % Line 8 has the units (TODO make this safer)
                units = strsplit(lineOfText,["\t"," "],'CollapseDelimiters',true);
                results(i).units = string(units(2:end-1));
            elseif lineNumber > 8 
                dataStr = strsplit(lineOfText,["\t"," "],'CollapseDelimiters',true);
                results(i).data(dataRow,:) = str2double(dataStr(2:end));
                dataRow = dataRow + 1;
            end
            lineOfText = fgetl(fid);
            lineNumber = lineNumber + 1;
        end
        fclose(fid);
    end
    cd(currDir);
end

function s = getSimulationCasesToRun(runname)
    switch runname
        case "compare2zeroPitch"
            sims(1).name = "still_p0";
            sims(2).name = "towSpeed05_p20"; 
            sims(3).name = "towSpeed1_p20"; 
            sims(4).name = "towSpeed15_p20"; 
            sims(5).name = "towSpeed2_p20"; 
            sims(6).name = "towSpeed25_p20"; 
            sims(7).name = "towSpeed3_p20"; 
            sims(8).name = "towSpeed35_p20";
        case "compare2twentyPitch"
            sims(1).name = "still_p0";
            sims(2).name = "towSpeed05_p20"; 
            sims(3).name = "towSpeed1_p20"; 
            sims(4).name = "towSpeed15_p20"; 
            sims(5).name = "towSpeed2_p20"; 
            sims(6).name = "towSpeed25_p20"; 
            sims(7).name = "towSpeed3_p20"; 
            sims(8).name = "towSpeed35_p20";
        otherwise
            error('unknown run name');
    end
end
% towedTurbineAerodyn.m
% Runs aerodyn analyses for the towed turbine project.
% Written by Rodney Metoyer
% All rights reserved


%% Clear the workspace
clearvars; close all; clc;
startingDir = pwd;

%% Set up simulations
runAerodynSimulations = true;
parseAerodynResults = true;
makeIndividualSimPlots = true;

% This set assumes a base (freestream) flow of 1.5 m/s. Uses a 1-meter
% radius (2m diametsr) blade with S814 airfoil. 675 simulations.
baseflow = 2.0;
maxTowSpeed = 10.0;
pitch = 0:10:0;
maxRelSpeed = baseflow+maxTowSpeed;
towspeeds = 8:0.5:maxRelSpeed-baseflow;

%simulationSet = "NCWT_standardTwist";
%simulationSet = "NCWT_lowTwist";
%simulationSet = "NCWT_veryLowTwist";
simulationSet = "NCWT_OppVeryLowTwist";
for i=1:1:numel(pitch)
    close all;
    for j=1:1:numel(towspeeds)
        towspeedsstr = string(towspeeds(j));
        towspeedsstr = strrep(towspeedsstr,".","_");
        sims(j).name = strcat(simulationSet,string(pitch(i)),"_towSpeed",towspeedsstr);
    end
    if ~strcmp(pwd,startingDir)
        error('Looks like a previous run was interupted. Wrong working folder.');
    end
    simGroupDirName = strcat(simulationSet, string(pitch(i)));
    runTheSimulationsGroup(sims,runAerodynSimulations,parseAerodynResults,baseflow,char(simGroupDirName),makeIndividualSimPlots);
    clear sims simGroupDirName
end


function runTheSimulationsGroup(sims,runAerodynSimulations,parseAerodynResults,flowspeed,simGroupDirName,makeIndividualSimPlots)
    %% Run AeroDyn
    if(runAerodynSimulations)
        currDir = cd([pwd '\simulationInputFiles']);
        for i=1:1:numel(sims)
            driverFileNames = [char(sims(i).name) '_driver.dvr'];
            command=['..\..\lib\AeroDyn_Driver_Win32 ' driverFileNames];
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
                    elseif strcmp('RtAeroMxh',sims(i).results(j).headers(k))
                        sims(i).RtAeroMxh(j) = mean(sims(i).results(j).data(:,k));
                    end
                    %sims(i).(sims(i).results(j).headers(k)) = sims(i).results(j).data(:,k);
                end
            end
            towspeed = unique(towspeed);
            if numel(towspeed)>1
                error('Should all be the same speed for every simulation case.');
            end
            sims(i).towspeed = towspeed;
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
    %syms a
    for i=1:1:numel(sims)
        %for j=1:1:numel(sims(i).RtAeroCp)
            %cp = sims(i).RtAeroCp(j);
            %S = vpasolve(4*a*(1-a)^2 == cp, a);
            %S = double(S);
            %S = S(S<=1.0);
            %if cp > 0
            %    S = S(S>=0.0);
            %end
            %if numel(S) > 1
            %    S = S(S<0.5);
            %    if numel(S)>1
            %        warning(['Multiple valid induction factor solutions. Check results.']);
            %        S
            %        S = S(1);                
            %    end
            %end
            %if numel(S) ~= 1
            %    disp('empty');
            %end
            %sims(i).anum(j) = S;
        %end
        % max power
        [sims(i).maxCp,ind] = max(sims(i).RtAeroCp);
        [sims(i).maxPower,ind2] = max(sims(i).RtAeroPwr);
        if ind ~= ind2
            error('This is a problem. They should be in the same place');
        end
        sims(i).thrustMaxCp = sims(i).RtAeroFxh(ind);
        %sims(i).indFactMaxCp = sims(i).anum(ind);

        casename{i} = char(sims(i).name);
        power(i) = sims(i).maxPower;
        thrust(i) = sims(i).thrustMaxCp;
        %induction(i) = sims(i).indFactMaxCp;
        towPower(i) = sims(i).thrustMaxCp*sims(i).towspeed;
        netPower(i) = sims(i).maxPower - sims(i).thrustMaxCp*sims(i).towspeed;
        % warning, the next line assumes still is the first simulation
        powerGain(i) = netPower(i) - sims(1).maxPower;
        towspeed(i) = sims(i).towspeed;

        disp([char(sims(i).name) ' Power: ' num2str(power(i),'%6.2f') ' (W)']);
        disp([char(sims(i).name) ' Thrust: ' num2str(thrust(i),'%6.2f') ' (N)']);
        %disp([char(sims(i).name) ' Induction: ' num2str(induction(i),'%6.4f') ' (unitless)']);
        disp([char(sims(i).name) ' Tow Power: ' num2str(towPower(i),'%6.4f') ' (W)']);
        disp([char(sims(i).name) ' Net Power: ' num2str(netPower(i),'%6.4f') ' (W)']);
        disp([char(sims(i).name) ' Power Gain: ' num2str(powerGain(i),'%6.4f') ' (W)']);
        disp([char(sims(i).name) ' Power Gain %: ' num2str(powerGain(i)/sims(1).maxPower*100,'%6.4f') ' (W)']);
        disp(' ');
    end
    
    power = power.'; thrust = thrust.'; %induction = induction.'; 
    towPower = towPower.';
    netPower = netPower.'; powerGain = powerGain.'; towspeed = towspeed.';
    %T = table(towspeed,power,thrust,induction,towPower,netPower,powerGain,'RowNames',casename);
    T = table(towspeed,power,thrust,towPower,netPower,powerGain,'RowNames',casename);

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
    
    if makeIndividualSimPlots
        for i=1:1:numel(sims)
            hfig = figure;
            ax = axes('Parent',hfig);
            plot(ax,sims(i).RtTSR,sims(i).RtAeroPwr,'*');
            xlabel(ax,'TSR'); ylabel(ax,'Power (W)');
            title(ax,sims(i).name,'Interpreter','none');
            ax.Color = 'none';
            export_fig(hfig,[writeDir '\figures\' char(sims(i).name) 'pwrVtsr.png'],'-transparent','-m3');

            hfig = figure;
            ax = axes('Parent',hfig);
            plot(ax,sims(i).RtTSR,sims(i).RtAeroCp,'*');
            xlabel(ax,'TSR'); ylabel(ax,'Power Coefficient (c_P)');
            title(ax,sims(i).name,'Interpreter','none');
            ax.Color = 'none';
            export_fig(hfig,[writeDir '\figures\' char(sims(i).name) 'CpVtsr.png'],'-transparent','-m3');

            hfig = figure;
            ax = axes('Parent',hfig);
            plot(ax,sims(i).RtTSR,sims(i).RtAeroFxh,'*');
            xlabel(ax,'TSR'); ylabel(ax,'Thrust (N)');
            title(ax,sims(i).name,'Interpreter','none');
            ax.Color = 'none';
            export_fig(hfig,[writeDir '\figures\' char(sims(i).name) 'ThrustVtsr.png'],'-transparent','-m3');

%             hfig = figure;
%             ax = axes('Parent',hfig);
%             plot(ax,sims(i).RtTSR,sims(i).anum,'*');
%             xlabel(ax,'TSR'); ylabel(ax,'Induction Factor');
%             title(ax,sims(i).name,'Interpreter','none');    
%             ax.Color = 'none';
%             export_fig(hfig,[writeDir '\figures\' char(sims(i).name) 'InductionVtsr.png'],'-transparent','-m3');
        end
    end
    
    disp(['Making Power Gain vs Tow Speed Figure for ' simGroupDirName]);
    hfig = figure;
    ax = axes('Parent',hfig);
    plot(ax,-towspeed,powerGain,'*');
    grid(ax,'on');
    xlabel(ax,'Tow Speed (m/s)'); ylabel(ax,'Power Gain (W)');
    title('Power Gain vs. Tow Speed','Interpreter','none');
    ax.Color = 'none';
    export_fig(hfig,[writeDir '\figures\powerGainVtowSpeed.png'],'-transparent','-m3');
    saveas(hfig,[writeDir '\figures\powerGainVtowSpeed.fig']);
    
    disp(['Making Net Power vs Tow Speed Figure for ' simGroupDirName]);
    hfig = figure;
    ax = axes('Parent',hfig);
    plot(ax,-towspeed,netPower,'*');
    grid(ax,'on');
    xlabel(ax,'Tow Speed (m/s)'); ylabel(ax,'Net Power (W)');
    %title('Net Power vs. Tow Speed','Interpreter','none');
    ax.Color = 'none';
    export_fig(hfig,[writeDir '\figures\netpowerVtowSpeed.png'],'-transparent','-m3');
    saveas(hfig,[writeDir '\figures\netpowerVtowSpeed.fig']);
end

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

%% copy-pasted stuff that I want to keep to save typing maybe
%sims(1).name = "Bahaj_MHK"; % Verification test that comes with Aerodyn
%sims(2).name = "Doman_MHK"; % Verification test that comes with Aerodyn
% sims(1).name = "still_p0";           sims(1).towspeed = 0;
% sims(2).name = "towSpeed3_p0";       sims(2).towspeed = 3;
% sims(3).name = "towSpeed3_p20";      sims(3).towspeed = 3;
% sims(4).name = "towSpeed6_p0";       sims(4).towspeed = 6;
% sims(5).name = "towSpeed6_p20";      sims(5).towspeed = 6;
% sims(6).name = "towSpeed3_p0_tsr2";  sims(6).towspeed = 3;
% sims(7).name = "towSpeed6_p0_tsr2";  sims(7).towspeed = 6;
% sims(8).name = "towSpeed3_p10_tsr2"; sims(8).towspeed = 3;
% sims(9).name = "towSpeed6_p10_tsr2"; sims(9).towspeed = 6;

% sims(1).name = "still_p0";           sims(1).towspeed = 0;
% sims(2).name = "towSpeed3_p20";      sims(2).towspeed = 3;

% sims(1).name = "towSpeed3_p20_tsr2";      sims(1).towspeed = 3;
% sims(2).name = "towSpeed6_p20_tsr2";      sims(2).towspeed = 6;

% sims(1).name = "towSpeed1_p10_tsr2"; sims(1).towspeed = 1;
% sims(2).name = "towSpeed2_p10_tsr2"; sims(2).towspeed = 2;
% sims(3).name = "towSpeed4_p10_tsr2"; sims(3).towspeed = 4;
% sims(4).name = "towSpeed5_p10_tsr2"; sims(4).towspeed = 5;

% sims(1).name = "towSpeed25_p10_tsr2"; sims(1).towspeed = 1;
% sims(2).name = "towSpeed35_p10_tsr2"; sims(2).towspeed = 2;

% simGroupName = 'Compare to bad stationary'; simGroupDirName = 'comp2bad';
% %simGroupName = 'Compare to good stationary'; simGroupDirName = 'comp2good';
% sims(1).name = "still_p20";
% %sims(1).name = "still_p0";
% sims(2).name = "towSpeed05_p20"; 
% sims(3).name = "towSpeed1_p20"; 
% sims(4).name = "towSpeed15_p20"; 
% sims(5).name = "towSpeed2_p20"; 
% sims(6).name = "towSpeed25_p20"; 
% sims(7).name = "towSpeed3_p20"; 
% sims(8).name = "towSpeed35_p20";

% simGroupName = 'Compare to good stationary no pitch'; simGroupDirName = 'compNoPitch';
% sims(1).name = "still_p0";
% sims(2).name = "towSpeed05_p0"; 
% sims(3).name = "towSpeed1_p0"; 
% sims(4).name = "towSpeed15_p0"; 
% sims(5).name = "towSpeed2_p0"; 
% sims(6).name = "towSpeed25_p0"; 
% sims(7).name = "towSpeed3_p0"; 
% sims(8).name = "towSpeed35_p0";
% sims(9).name = "towSpeed4_p0";

% simGroupName = 'Compare to good stationary no pitch'; simGroupDirName = 'compGood2tsr';
% sims(1).name = "still_p0";
% sims(2).name = "towSpeed05_p20_2tsr"; 
% sims(3).name = "towSpeed1_p20_2tsr"; 
% sims(4).name = "towSpeed15_p20_2tsr"; 
% sims(5).name = "towSpeed2_p20_2tsr"; 
% sims(6).name = "towSpeed25_p20_2tsr"; 
% sims(7).name = "towSpeed3_p20_2tsr"; 
% sims(8).name = "towSpeed35_p20_2tsr";
% sims(9).name = "towSpeed4_p20_2tsr";

% simGroupName = 'Compare to good stationary no pitch'; simGroupDirName = 'compGood2tsrP15';
% sims(1).name = "still_p0";
% sims(2).name = "towSpeed05_p15_2tsr"; 
% sims(3).name = "towSpeed1_p15_2tsr"; 
% sims(4).name = "towSpeed15_p15_2tsr"; 
% sims(5).name = "towSpeed2_p15_2tsr"; 
% sims(6).name = "towSpeed25_p15_2tsr"; 
% sims(7).name = "towSpeed3_p15_2tsr"; 
% sims(8).name = "towSpeed35_p15_2tsr";
% sims(9).name = "towSpeed4_p15_2tsr";

% simGroupName = 'Compare to good stationary no pitch'; simGroupDirName = 'compGood2tsrP10';
% sims(1).name = "still_p0";
% sims(2).name = "towSpeed05_p10_2tsr"; 
% sims(3).name = "towSpeed1_p10_2tsr"; 
% sims(4).name = "towSpeed15_p10_2tsr"; 
% sims(5).name = "towSpeed2_p10_2tsr"; 
% sims(6).name = "towSpeed25_p10_2tsr"; 
% sims(7).name = "towSpeed3_p10_2tsr"; 
% sims(8).name = "towSpeed35_p10_2tsr";
% sims(9).name = "towSpeed4_p10_2tsr";

% simGroupName = 'Compare to good stationary no pitch'; simGroupDirName = 'compGood2tsrP5';
% sims(1).name = "still_p0";
% sims(2).name = "towSpeed05_p5_2tsr"; 
% sims(3).name = "towSpeed1_p5_2tsr"; 
% sims(4).name = "towSpeed15_p5_2tsr"; 
% sims(5).name = "towSpeed2_p5_2tsr"; 
% sims(6).name = "towSpeed25_p5_2tsr"; 
% sims(7).name = "towSpeed3_p5_2tsr"; 
% sims(8).name = "towSpeed35_p5_2tsr";
% sims(9).name = "towSpeed4_p5_2tsr";

% simGroupName = 'Compare to good stationary no pitch'; simGroupDirName = 'compGood2tsrP0';
% sims(1).name = "still_p0";
% sims(2).name = "towSpeed05_p0_2tsr"; 
% sims(3).name = "towSpeed1_p0_2tsr"; 
% sims(4).name = "towSpeed15_p0_2tsr"; 
% sims(5).name = "towSpeed2_p0_2tsr"; 
% sims(6).name = "towSpeed25_p0_2tsr"; 
% sims(7).name = "towSpeed3_p0_2tsr"; 
% sims(8).name = "towSpeed35_p0_2tsr";
% sims(9).name = "towSpeed4_p0_2tsr";

% simGroupName = 'Compare to good stationary no pitch'; simGroupDirName = 'compGood1tsrP0';
% sims(1).name = "still_p0";
% sims(2).name = "towSpeed05_p0_1tsr"; 
% sims(3).name = "towSpeed1_p0_1tsr"; 
% sims(4).name = "towSpeed15_p0_1tsr"; 
% sims(5).name = "towSpeed2_p0_1tsr"; 
% sims(6).name = "towSpeed25_p0_1tsr"; 
% sims(7).name = "towSpeed3_p0_1tsr"; 
% sims(8).name = "towSpeed35_p0_1tsr";
% sims(9).name = "towSpeed4_p0_1tsr";

% simGroupName = 'Compare to good stationary no pitch'; simGroupDirName = 'compGood1tsrP5';
% sims(1).name = "still_p0";
% sims(2).name = "towSpeed05_p5_1tsr"; 
% sims(3).name = "towSpeed1_p5_1tsr"; 
% sims(4).name = "towSpeed15_p5_1tsr"; 
% sims(5).name = "towSpeed2_p5_1tsr"; 
% sims(6).name = "towSpeed25_p5_1tsr"; 
% sims(7).name = "towSpeed3_p5_1tsr"; 
% sims(8).name = "towSpeed35_p5_1tsr";
% sims(9).name = "towSpeed4_p5_1tsr";

% pitchstr = "10";
% simGroupDirName = 'compGood1tsrP10';
% speedsstr = ["0","05","1","15","2","25","3","35","4"];
% sims(1).name = "still_p0";
% for i = 2:1:9
%     sims(i).name = strcat("towSpeed",speedsstr(i),"_p",pitchstr,"_1tsr");
% end
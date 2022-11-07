% makeAerodynInputFiles.m
% script to make aerodyn input files
clearvars; close all; clc;


pitch = 0;         % Array of pitch values to try
baseFlow = 2.0;    % This is the freestream flow speed in m/s
topTowSpeed = 20;  % Maximum tow speed 
stopTSR = 3.0;     % Highest TSR that you want

flowspeed = baseFlow:0.5:topTowSpeed+baseFlow; % This is the current speed plus the tow speed
hubRadius = 0.01481667; % This goes with the 35.56 cm diameter blade defined by ncsuwt_blade.dat
%bladeFileName = 'ncsuwt_blade_lowTwist.dat';       %%%%% REMEMBER TO CHANGE THE HUB RADIUS IF YOU USE A DIFFERENT ROTOR %%%%%%
bladeFileName = 'ncsuwt_blade_OppVeryLowTwist.dat'; %%%%% REMEMBER TO CHANGE THE HUB RADIUS IF YOU USE A DIFFERENT ROTOR %%%%%%
%bladeFileName = 'ncsuwt_blade_OppVeryLowTwist.dat';
runname = "NCWT_OppVeryLowTwist";
%airfoil = '"Airfoils\SG6040.dat"';
airfoil = '"Airfoils\NRELS814.dat"';

% WARNING - Still need to change fluid parameters!!!!!

for i=1:1:numel(pitch)
    for j=1:1:numel(flowspeed)
        towspeedsstr = string(flowspeed(j)-baseFlow);
        towspeedsstr = strrep(towspeedsstr,".","_");
        %p.runName = strcat("MB_F1d5_SG6040_TSR2","_p",string(pitch(i)),"_towSpeed",towspeedsstr);
        %p.runName = strcat("MB_F1d5_SG6040_p",string(pitch(i)),"_towSpeed",towspeedsstr);
        %p.runName = strcat("EastRiver_2_0_p",string(pitch(i)),"_towSpeed",towspeedsstr);
        %p.runName = strcat("NCWT_standardTwist",string(pitch(i)),"_towSpeed",towspeedsstr);
        %p.runName = strcat("NCWT_veryLowTwist",string(pitch(i)),"_towSpeed",towspeedsstr);
        p.runName = strcat(runname,string(pitch(i)),"_towSpeed",towspeedsstr);
        p.pitch = pitch(i);
        p.flowspeed = flowspeed(j);
        p.stopTSR = stopTSR;
        p.hubRadius = hubRadius;
        makeAerodynInputFilesSF(p,bladeFileName,airfoil);
    end
end



function makeAerodynInputFilesSF(inStruct,bfn,afn)
    % INPUT: inParam - a data structure to use for inputs
    if nargin < 3
        bfn = 'meter_blade.dat';
        afn = '"Airfoils\NRELS814.dat"';
    end

    % runName is the name of the simulation case (the "run"). You will get
    % two aerodyn input files out of this. ONe called <runName>_input.dat
    % and one named <runName>_driver.dvr
    
    pitch = inStruct.pitch;
    runName = char(inStruct.runName);
    flowspeed = inStruct.flowspeed;
    %bladeFileName = 'meter_blade_tsr1.dat';
    %bladeFileName = 'meter_blade_tsr2.dat';
    %bladeFileName = 'meter_blade.dat';
    %bladeFileName = 'ncsuwt_blade.dat';
    bladeFileName = bfn;
    TSR = linspace(0,inStruct.stopTSR,51);
    %rotspeed = TSR*60*flowspeed/(2*pi*1);
    rotspeed = TSR*60*flowspeed/(2*pi*0.1778); % Last number is radius. This is for the experimental rotor defined by ncsuwt_blade.dat

    % Setup driver stuff
    dp.runName = runName;
    % driverParams.numCases = 32;
    % Number of 
    dp.numBlades = 3;     %1,2,or 3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %dp.hubRad   = 0.0833;   %Hub radius in meters
    dp.hubRad   = inStruct.hubRadius; 
    dp.hubHt    = 100.0;    %Hub height in meters (off the ground for ground interaction)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dp.overhang = 0.0;
    dp.shftTilt = 0.0;
    dp.precone  = 0.0;

    yawindx = 0;
    numcases = 0;
    for yaw = 0:1:0.1
        yawindx = yawindx+1;
        rotindx = 0;
        for i=1:1:numel(rotspeed)
            rotindx = rotindx + 1;
            % todo make a flowspeed loop as well
            numcases = numcases + 1;

            dp.flowSpeed(numcases) = flowspeed;
            dp.shearExp(numcases) = 0.0;
            dp.rotSpeed(numcases) = rotspeed(i); % rpm
            dp.pitch(numcases) = pitch; %-6;    %Pitch in degrees
            dp.yaw(numcases) = yaw;           %Yaw in degrees
            dp.dt(numcases) = 0.1;            %Timestep size in seconds
            dp.tmax(numcases) = 5.6;         %Max runtime in seconds 
        end
    end
    dp.numCases = numcases;
    dp.totalCases = numcases;
    dp.yawCases = yawindx;
    dp.TSRCases = rotindx;

    % Make the blade file

    %bladeFileName = [char(runNames(1)) '_blade.dat'];
    %MakeBladeFile(bladeFileName)

    % Make the input file
    ip.inputFileName = [pwd '\simulationInputFiles\' dp.runName '_input.dat'];
    %inptFileName = ip.inputFileName; %obsolete
    ip.bladeFileName = bladeFileName; % TODO make string array - input file supports up to 3 different blades with different blade files.
    ip.echo = 'False';          % Echo the input to "<rootname>.AD.ech"?
    ip.WakeMod = 1;             % Type of wake/induction model (switch) {0=none, 1=BEMT}
    ip.AFAeroMod = 1;           % Type of blade airfoil aerodynamics model (switch) {1=steady model, 2=Beddoes-Leishman unsteady model} [must be 1 when linearizing]
    ip.CavitCheck = 'False';    % Perform cavitation check?
    %ip.AirDens = 998.0;        % Fluiddensity (kg/m^3)
    ip.AirDens = 1.225;        % Fluiddensity (kg/m^3)
    %ip.KinVisc = 1.0E-3;        % Kinematic viscosity (m^2/s)
    ip.KinVisc = 1.48E-5;        % Kinematic viscosity (m^2/s)
    ip.SpdSound = 1500.0;       % Speed of sound (m/s)
    ip.Patm = 103500.0;         % Atmospheric pressure (Pa) [used only when CavitCheck=True]
    ip.Pvap = 1700.0;           % Vapour pressure of fluid (Pa) [used only when CavitCheck=True]
    ip.FluidDepth = 0.1;        % Water depth above mid-hub height (m) [used only when CavitCheck=True]
    ip.SkewMod = 1;             % Type of skewed-wake correction model (switch) {1=uncoupled, 2=Pitt/Peters, 3=coupled} [used only when WakeMod=1]
    ip.TipLoss = 'False';       % Use the Prandtl tip-loss model? (flag) [used only when WakeMod=1]
    ip.HubLoss = 'False';       % Use the Prandtl hub-loss model? (flag) [used only when WakeMod=1]
    ip.TanInd = 'False';        % Include tangential induction in BEMT calculations? (flag) [used only when WakeMod=1]
    ip.AIDrag = 'False';        % Include the drag term in the axial-induction calculation? (flag) [used only when WakeMod=1]
    ip.TIDrag = 'False';        % Include the drag term in the tangential-induction calculation? (flag) [used only when WakeMod=1 and TanInd=TRUE]
    %ip.airfoil = '"Airfoils\SG6040.dat"';
    %ip.airfoil = '"Airfoils\NRELS814.dat"'; 
    %ip.airfoil = '"Airfoils\DU30_A17_100Hz.dat"';
    ip.airfoil = afn;
    % Need to specifiy columns in the airfoil file for coeffs. 0 = not there
    ip.InCol_Alfa = 1;          % The column in the airfoil tables that contains the angle of attack
    ip.InCol_Cl = 2;            % The column in the airfoil tables that contains the lift coefficient
    ip.InCol_Cd = 3;            % The column in the airfoil tables that contains the drag coefficient
    ip.InCol_Cm = 0;            % The column in the airfoil tables that contains the pitching-moment coefficient; use zero if there is no Cm column
    ip.InCol_Cpmin = 0;         % The column in the airfoil tables that contains the Cpmin coefficient; use zero if there is no Cpmin column
    ip.NumAFfiles = 1;          % Number of airfoil files used
    ip.outputs = ["RtSpeed","RtTSR","RtAeroPwr","RtAeroFxh","RtAeroCp",...
        "RtAeroCt","B1N7Alpha","B1N7Cl","B1N7Cd","B1N7Cm",...
        "B1Azimuth","B2Azimuth","B1Pitch","RtAeroFyh","RtAeroFzh",...
        "RtAeroMxh","RtAeroMyh","RtAeroMzh"]; %TODO this can be a gui by iteself maybe
    makeInputFile(ip);

    % Make the driver file
    % Name should be runName_driver
    dp.bladeFileName = bladeFileName;
    dp.inputFileName = ip.inputFileName;
    dp.driverFileName = [pwd '\simulationInputFiles\' dp.runName '_driver.dvr'];
    %runName = driverParams.driverFileName; %obsolete
    makeDriverFile(dp);
end

function [] = makeDriverFile(dp)
    fid = fopen(dp.driverFileName,'w');
    fprintf(fid,'%s\r\n','-----  AeroDyn Driver v15.04.02a Input File  -----------------------------------');
    fprintf(fid,'%s\r\n',dp.bladeFileName);
    fprintf(fid,'%s\r\n','-----  Input Configuration  ----------------------------------------------------');
    fprintf(fid,'%s\r\n','FALSE          Echo');
    fprintf(fid,'%s\r\n',['"' dp.inputFileName '"']);
    fprintf(fid,'%s\r\n','-----  Turbine Data  -----------------------------------------------------------');
    fprintf(fid,'%s %d %s\r\n','        ',dp.numBlades, 'NumBlades');
    fprintf(fid,'%s %6.4f %s\r\n','        ',dp.hubRad, 'HubRad');
    fprintf(fid,'%s %3.2f %s\r\n','        ',dp.hubHt, 'HubHt');
    fprintf(fid,'%s %3.2f %s\r\n','        ',dp.overhang, 'Overhang');
    fprintf(fid,'%s %3.2f %s\r\n','        ',dp.shftTilt, 'ShftTilt');
    fprintf(fid,'%s %3.2f %s\r\n','        ',dp.precone, 'Precone');
    fprintf(fid,'%s\r\n','-----  I/O Settings  -----------------------------------------------------------');
    fprintf(fid,'%s\r\n',[dp.runName '            OutFileRoot']);
    fprintf(fid,'%s\r\n','True          TabDel');
    fprintf(fid,'%s\r\n','"ES15.6E3"    OutFmt');
    fprintf(fid,'%s\r\n','True          Beep');
    fprintf(fid,'%s\r\n','-----  Combined-Case Analysis  -------------------------------------------------');
    fprintf(fid,'%s\r\n',[num2str(dp.numCases) '    NumCases']);
    fprintf(fid,'%s\r\n','FlowSpeed   ShearExp   RotSpd   Pitch   Yaw   dT   Tmax');
    fprintf(fid,'%s\r\n','(m/s)       (-)        (rpm)    (deg)   (deg) (s)  (s)');
    %fprintf(fid,'%.7e %.7e %4.1f %.7e %.7e %4.3f %4.3f\r\n',dp.flowSpeed,dp.shearExp,dp.rotSpeed,dp.pitch,dp.yaw,dp.dt,dp.tmax);
    for i=1:1:dp.numCases
        fprintf(fid,'%.7e %.7e %4.1f %.7e %.7e %4.3f %4.3f\r\n',dp.flowSpeed(i),dp.shearExp(i),dp.rotSpeed(i),dp.pitch(i),dp.yaw(i),dp.dt(i),dp.tmax(i));
    end
    fid = fclose(fid);
end

function [] = makeInputFile(ip)
    try
        fid = fopen(ip.inputFileName,'w');
    catch e
        disp(e.message);
    end
    fprintf(fid,'%s\r\n','-----  AeroDyn v15.04.02a Input File  -----------------------------------');
    fprintf(fid,'%s\r\n',ip.bladeFileName);
    fprintf(fid,'%s\r\n','======  General Options  ============================================================================');
    fprintf(fid,'%s %s\r\n',ip.echo,'        Echo');
    fprintf(fid,'%s\r\n','"default"     DTAero');
    fprintf(fid,'%s %d %s\r\n','        ',ip.WakeMod,'   WakeMod');
    fprintf(fid,'%s %d %s\r\n','        ',ip.AFAeroMod,'   AFAeroMod');
    fprintf(fid,'%s\r\n','         0    TwrPotent');
    fprintf(fid,'%s\r\n','False         TwrShadow');
    fprintf(fid,'%s\r\n','False         TwrAero');
    fprintf(fid,'%s\r\n','False         FrozenWake');
    fprintf(fid,'%s %s\r\n',ip.CavitCheck,'         CavitCheck');
    fprintf(fid,'%s\r\n','======  Environmental Conditions  ===================================================================');
    fprintf(fid,'%s %d %s\r\n','          ',ip.AirDens,'           AirDens');
    fprintf(fid,'%s %d %s\r\n','          ',ip.KinVisc,'   KinVisc');
    fprintf(fid,'%s %d %s\r\n','          ',ip.SpdSound,'           SpdSound');
    fprintf(fid,'%s %d %s\r\n','          ',ip.Patm,'         Patm');
    fprintf(fid,'%s %d %s\r\n','          ',ip.Pvap,'           Pvap');
    fprintf(fid,'%s %d %s\r\n','          ',ip.FluidDepth,'   FluidDepth');
    fprintf(fid,'%s\r\n','======  Blade-Element/Momentum Theory Options  ====================================================== [used only when WakeMod=1]');
    fprintf(fid,'%d %s\r\n',ip.SkewMod,'   SkewMod');
    fprintf(fid,'%s %s\r\n',ip.TipLoss,'         TipLoss');
    fprintf(fid,'%s %s\r\n',ip.HubLoss,'         HubLoss');
    fprintf(fid,'%s %s\r\n',ip.TanInd,'         TanInd');
    fprintf(fid,'%s %s\r\n',ip.AIDrag,'         AIDrag');
    fprintf(fid,'%s %s\r\n',ip.TIDrag,'         TIDrag');
    fprintf(fid,'%s\r\n','"default"     IndToler');
    fprintf(fid,'%s\r\n','10000   MaxIter');
    fprintf(fid,'%s\r\n','======  Beddoes-Leishman Unsteady Airfoil Aerodynamics Options  ===================================== [used only when AFAeroMod=2]');
    fprintf(fid,'%s\r\n','          2   UAMod');
    fprintf(fid,'%s\r\n','True          FLookup');
    fprintf(fid,'%s\r\n','======  Airfoil Information =========================================================================');
    fprintf(fid,'%d %s\r\n',ip.InCol_Alfa,'   InCol_Alfa         - The column in the airfoil tables that contains the angle of attack (-)');
    fprintf(fid,'%d %s\r\n',ip.InCol_Cl,'   InCol_Cl           - The column in the airfoil tables that contains the lift coefficient (-)');
    fprintf(fid,'%d %s\r\n',ip.InCol_Cd,'   InCol_Cd           - The column in the airfoil tables that contains the drag coefficient (-)');
    fprintf(fid,'%d %s\r\n',ip.InCol_Cm,'   InCol_Cm           - The column in the airfoil tables that contains the pitching-moment coefficient; use zero if there is no Cm column (-)');
    fprintf(fid,'%d %s\r\n',ip.InCol_Cpmin,'   InCol_Cpmin        - The column in the airfoil tables that contains the Cpmin coefficient; use zero if there is no Cpmin column (-)');
    fprintf(fid,'%d %s\r\n',ip.NumAFfiles,'   NumAFfiles         - Number of airfoil files used (-)');
    fprintf(fid,'%s %s\r\n',ip.airfoil,'    AFNames            - Airfoil file names (NumAFfiles lines) (quoted strings)');
    fprintf(fid,'%s\r\n','======  Rotor/Blade Properties  =====================================================================');
    fprintf(fid,'%s\r\n','True          UseBlCm            - Include aerodynamic pitching moment in calculations?  (flag)');
    fprintf(fid,'%s%s%s %s\r\n','"',ip.bladeFileName,'"','                 ADBlFile(1)        - Name of file containing distributed aerodynamic properties for Blade #1 (-)');
    fprintf(fid,'%s%s%s %s\r\n','"',ip.bladeFileName,'"','                 ADBlFile(2)        - Name of file containing distributed aerodynamic properties for Blade #2 (-) [unused if NumBl < 2]');
    fprintf(fid,'%s%s%s %s\r\n','"',ip.bladeFileName,'"','                 ADBlFile(3)        - Name of file containing distributed aerodynamic properties for Blade #3 (-) [unused if NumBl < 3]');
    fprintf(fid,'%s\r\n','======  Tower Influence and Aerodynamics ============================================================= [used only when TwrPotent/=0, TwrShadow=True, or TwrAero=True]');
    fprintf(fid,'%s\r\n','         1    NumTwrNds         - Number of tower nodes used in the analysis  (-) [used only when TwrPotent/=0, TwrShadow=True, or TwrAero=True]');
    fprintf(fid,'%s\r\n','TwrElev        TwrDiam        TwrCd');
    fprintf(fid,'%s\r\n','(m)              (m)           (-)');
    fprintf(fid,'%s\r\n','1.0000000E+00  6.0000000E+00  1.0000000E+00  ');
    fprintf(fid,'%s\r\n','======  Outputs  ====================================================================================');
    fprintf(fid,'%s\r\n','True         SumPrint            - Generate a summary file listing input options and interpolated properties to "<rootname>.AD.sum"?  (flag)');
    fprintf(fid,'%s\r\n','          9   NBlOuts             - Number of blade node outputs [0 - 9] (-)');
    fprintf(fid,'%s\r\n','          1,	2,	3,	4,          5,          6,   	7,       8,          9    BlOutNd             - Blade nodes whose values will be output  (-)');
    fprintf(fid,'%s\r\n','          0   NTwOuts             - Number of tower node outputs [0 - 9]  (-)');
    fprintf(fid,'%s\r\n','          1,          2,          3,          4,          5    TwOutNd             - Tower nodes whose values will be output  (-)');
    fprintf(fid,'%s\r\n','OutList             - The next line(s) contains a list of output parameters.  See OutListParameters.xlsx for a listing of available output channels, (-)');
    for i=1:1:length(ip.outputs)
        fprintf(fid,'%s%s%s\r\n','"',ip.outputs(i),'"');
    end
    fprintf(fid,'%s\r\n','END of input file (the word "END" must appear in the first 3 columns of this last OutList line)');
    fid = fclose(fid);
end
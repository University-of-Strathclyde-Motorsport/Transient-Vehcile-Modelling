%Read Vehicle Parameters Function
% - Murray
%This function creates an array of base parameters, no calculations or model generation is done here.
%The car 'model' is created in generatevehicle.
function [Vehicles] = readvehiclesheet(filename)

%% Read vehicle file
parameters = table2cell(read_parameters(filename,'Parameters'));
hpTorque = read_torque(filename,'High Power Map');
lpTorque = read_torque(filename,'Low Power Map');
swpits = str2double(string(parameters(3,2)));

%Create new arrays for each column. Not nice but matlab 
%wasn't playing well when working directly with parameters().
fieldnames = string(parameters(:,1));
mincol = str2double(string(parameters(:,2)));
maxcol = string(parameters(:,3));
maxcol = fillmissing(maxcol,'constant',"0");
maxcol = str2double(maxcol);
stringfields = ["Name","Version","FrontTyreName","RearTyreName","FuelType","DriveType"];

%Check there are only two Max values.
%Otherwise a 4d output would be required.
nmax = nnz(maxcol);
if nmax > 2
    error('Check input excel sheet. \n%s',...
        'There are currently more than 2 max fields.');
end

%Detect if user has asked for sweep but has not set limits or vice versa
if nmax == 0 && swpits > 1
    error("Sweep set but no field has a max value")
elseif swpits == 1 && nmax > 0
    error("Field limits set, but only one sweep iteration.")
end

%Define Vehicle array size.
rows = swpits;
if nmax == 2
    cols = swpits;
else
    cols = 1;
end

%Preallocate Vehicle array
Vehicles = cell(rows,cols);




%Create single structure containing arrays of sweeps for each field.
for r = 1:rows
    for c = 1:cols
        Vehicle = writerows(r,c);
        Vehicle.EngineSpeed = table2array(hpTorque(:,1));
        Vehicle.HighPowerTorque = table2array(hpTorque(:,2));
        Vehicle.LowPowerTorque = table2array(lpTorque(:,2));
        VP{r,c} = Vehicle;
    end
end

%Preallocate array of models
Vehicles = cell(size(VP,1),size(VP,2));

for r = 1:size(VP,1)
    for c = 1:size(VP,2)
        %Struct of current parameters
        vp = VP{r,c};

        %Name&Version
        % info
        veh.name = vp.Name;
        veh.type = vp.Type;
        if isfield(vp,'SweptParameter1') == 1
            veh.sweep1 = vp.SweptParameter1;
            veh.sweep1data = vp.SweptParameter1Data;
        end
        if isfield(vp,'SweptParameter2') == 1
            veh.sweep2 = vp.SweptParameter2;
            veh.sweep2data = vp.SweptParameter2Data;
        end

        % mass
        veh.M = vp.TotalMass;
        veh.df = vp.FrontMassDistribution/100;
        % wheelbase
        veh.L = vp.Wheelbase/1000;
        % steering rack ratio
        veh.rack = vp.SteeringRackRatio;
        % aerodynamics
        veh.Cl = vp.LiftCoefficientCL;
        veh.Cd = vp.DragCoefficientCD;
        veh.factor_Cl = vp.CLScaleMultiplier;
        veh.factor_Cd = vp.CDScaleMultiplier;
        veh.da = vp.FrontAeroDistribution/100;
        veh.A = vp.FrontalArea;
        veh.rho = vp.AirDensity;
        % brakes
        veh.br_disc_d = vp.DiscOuterDiameter/1000;
        veh.br_pad_h = vp.PadHeight/1000;
        veh.br_pad_mu = vp.PadFrictionCoefficient;
        veh.br_nop = vp.CaliperNumberofPistons;
        veh.br_pist_d = vp.CaliperPistonDiameter/1000;
        veh.br_mast_d = vp.MasterCylinderPistonDiameter/1000;
        veh.br_ped_r = vp.PedalRatio;
        % tyres
        veh.factor_grip = vp.GripFactorMultiplier;
        veh.tyre_radius = vp.TyreRadius/1000;
        veh.Cr = vp.RollingResistance;
        veh.mu_x = vp.LongitudinalFrictionCoefficient;
        veh.mu_x_M = vp.LongitudinalFrictionLoadRating;
        veh.sens_x = vp.LongitudinalFrictionSensitivity;
        veh.mu_y = vp.LateralFrictionCoefficient;
        veh.mu_y_M = vp.LateralFrictionLoadRating;
        veh.sens_y = vp.LateralFrictionSensitivity;
        veh.CF = vp.FrontCorneringStiffness;
        veh.CR = vp.RearCorneringStiffness;
        % engine
        veh.factor_power = vp.PowerFactorMultiplier;
        veh.n_thermal = vp.ThermalEfficiency;
        veh.fuel_LHV = vp.FuelLowerHeatingValue;
        % drivetrain
        veh.drive = vp.DriveType;
        veh.shift_time = vp.GearShiftTime;
        veh.n_primary = vp.PrimaryGearEfficiency/100;
        veh.n_final = vp.FinalGearEfficiency/100;
        veh.n_gearbox = vp.GearboxEfficiency/100;
        veh.ratio_primary = vp.PrimaryGearReduction;
        veh.ratio_final = vp.FinalGearReduction;
        veh.ratio_gearbox(1) = vp.Gear1Ratio;
        veh.ratio_gearbox(2) = vp.Gear2Ratio;
        veh.ratio_gearbox(3) = vp.Gear3Ratio;
        veh.ratio_gearbox(4) = vp.Gear4Ratio;
        veh.ratio_gearbox(5) = vp.Gear5Ratio;
        veh.ratio_gearbox(6) = vp.Gear6Ratio;
        veh.ratio_gearbox(7) = vp.Gear7Ratio;
        veh.ratio_gearbox(8) = vp.Gear8Ratio;
        veh.ratio_gearbox(9) = vp.Gear9Ratio;
        veh.ratio_gearbox(10) = vp.Gear10Ratio;
        for g = 1:length(veh.ratio_gearbox)
            if veh.ratio_gearbox(g) == 0
                veh.nog = g-1;
                veh.ratio_gearbox = veh.ratio_gearbox(veh.ratio_gearbox ~= 0);
                break
            end
        end
        veh.nog = length(veh.ratio_gearbox) ;

        %% Brake Model

        veh.br_pist_a = veh.br_nop*pi*(veh.br_pist_d/1000)^2/4 ; % [m2]
        veh.br_mast_a = pi*(veh.br_mast_d/1000)^2/4 ; % [m2]
        veh.beta = veh.tyre_radius/(veh.br_disc_d/2-veh.br_pad_h/2)/veh.br_pist_a/veh.br_pad_mu/4 ; % [Pa/N] per wheel
        veh.phi = veh.br_mast_a/veh.br_ped_r*2 ; % [-] for both systems


        %% Steering Model

        veh.a = (1-veh.df)*veh.L ; % distance of front axle from center of mass [mm]
        veh.b = -veh.df*veh.L ; % distance of rear axle from center of mass [mm]
        veh.C = 2*[veh.CF,veh.CF+veh.CR;veh.CF*veh.a,veh.CF*veh.a+veh.CR*veh.b] ; % steering model matrix

        %% Driveline Model

        % fetching engine curves
        veh.en_speed_curve = vp.EngineSpeed ; % [rpm]
        veh.en_torque_curve = vp.EngineTorque ; % [N*m]
        veh.en_power_curve = veh.en_torque_curve.*veh.en_speed_curve*2*pi/60 ; % [W]
        % memory preallocation
        % wheel speed per gear for every engine speed value
        veh.wheel_speed_gear = zeros(length(veh.en_speed_curve),veh.nog) ;
        % vehicle speed per gear for every engine speed value
        veh.vehicle_speed_gear = zeros(length(veh.en_speed_curve),veh.nog) ;
        % wheel torque per gear for every engine speed value
        veh.wheel_torque_gear = zeros(length(veh.en_torque_curve),veh.nog) ;
        % calculating values for each gear and engine speed
        for i=1:veh.nog
            veh.wheel_speed_gear(:,i) = veh.en_speed_curve/veh.ratio_primary/veh.ratio_gearbox(i)/veh.ratio_final ;
            veh.vehicle_speed_gear(:,i) = veh.wheel_speed_gear(:,i)*2*pi/60*veh.tyre_radius ;
            veh.wheel_torque_gear(:,i) = veh.en_torque_curve*veh.ratio_primary*veh.ratio_gearbox(i)*veh.ratio_final*veh.n_primary*veh.n_gearbox*veh.n_final ;
        end

        % minimum and maximum vehicle speeds
        veh.v_min = min(veh.vehicle_speed_gear,[],'all') ;
        veh.v_max = max(veh.vehicle_speed_gear,[],'all') ;

        % new speed vector for fine meshing
        dv = 0.5/3.6 ;
        veh.vehicle_speed = linspace(veh.v_min,veh.v_max,(veh.v_max-veh.v_min)/dv)' ;
        % memory preallocation
        % gear
        veh.gear = zeros(length(veh.vehicle_speed),1) ;
        % engine tractive force
        veh.fx_engine = zeros(length(veh.vehicle_speed),1) ;
        % engine tractive force per gear
        veh.fx = zeros(length(veh.vehicle_speed),veh.nog) ;
        % optimising gear selection and calculating tractive force
        for i=1:length(veh.vehicle_speed)
            % going through the gears
            for j=1:veh.nog
                veh.fx(i,j) = interp1(veh.vehicle_speed_gear(:,j),veh.wheel_torque_gear(:,j)/veh.tyre_radius,veh.vehicle_speed(i),'linear',0) ;
            end
            % getting maximum tractive force and gear
            [veh.fx_engine(i),veh.gear(i)] = max(veh.fx(i,:)) ;
        end
        % adding values for 0 speed to vectors for interpolation purposes at low speeds
        veh.vehicle_speed = [0;veh.vehicle_speed] ;
        veh.gear = [veh.gear(1);veh.gear] ;
        veh.fx_engine = [veh.fx_engine(1);veh.fx_engine] ;
        % final vectors
        % engine speed
        veh.engine_speed = veh.ratio_final*veh.ratio_gearbox(veh.gear)*veh.ratio_primary.*veh.vehicle_speed/veh.tyre_radius*60/2/pi ;
        % wheel torque
        veh.wheel_torque = veh.fx_engine*veh.tyre_radius ;
        % engine torque
        veh.engine_torque = veh.wheel_torque/veh.ratio_final./veh.ratio_gearbox(veh.gear)/veh.ratio_primary/veh.n_primary/veh.n_gearbox/veh.n_final ;
        % engine power
        veh.engine_power = veh.engine_torque.*veh.engine_speed*2*pi/60 ;

        %% Shifting Points and Rev Drops

        % finding gear changes
        gear_change = diff(veh.gear) ; % gear change will appear as 1
        % getting speed right before and after gear change
        gear_change = logical([gear_change;0]+[0;gear_change]) ;
        % getting engine speed at gear change
        engine_speed_gear_change = veh.engine_speed(gear_change) ;
        % getting shift points
        shift_points = engine_speed_gear_change(1:2:length(engine_speed_gear_change)) ;
        % getting arrive points
        arrive_points = engine_speed_gear_change(2:2:length(engine_speed_gear_change)) ;
        % calculating revdrops
        rev_drops = shift_points-arrive_points ;
        % creating shifting table
        rownames = cell(veh.nog-1,1) ;
        for i=1:veh.nog-1
            rownames(i) = {[num2str(i,'%d'),'-',num2str(i+1,'%d')]} ;
        end
        veh.shifting = table(shift_points,arrive_points,rev_drops,'RowNames',rownames) ;


        %% Force model

        % gravitational constant
        g = 9.81 ;
        % drive and aero factors
        switch veh.drive
            case 'RWD'
                veh.factor_drive = (1-veh.df) ; % weight distribution
                veh.factor_aero = (1-veh.da) ; % aero distribution
                veh.driven_wheels = 2 ; % number of driven wheels
            case 'FWD'
                veh.factor_drive = veh.df ;
                veh.factor_aero = veh.da ;
                veh.driven_wheels = 2 ;
            otherwise % AWD
                veh.factor_drive = 1 ;
                veh.factor_aero = 1 ;
                veh.driven_wheels = 4 ;
        end
        % Z axis
        veh.fz_mass = -veh.M*g ;
        veh.fz_aero = 1/2*veh.rho*veh.factor_Cl*veh.Cl*veh.A*veh.vehicle_speed.^2 ;
        veh.fz_total = veh.fz_mass+veh.fz_aero ;
        veh.fz_tyre = (veh.factor_drive*veh.fz_mass+veh.factor_aero*veh.fz_aero)/veh.driven_wheels ;
        % x axis
        veh.fx_aero = 1/2*veh.rho*veh.factor_Cd*veh.Cd*veh.A*veh.vehicle_speed.^2 ;
        veh.fx_roll = veh.Cr*abs(veh.fz_total) ;
        veh.fx_tyre = veh.driven_wheels*(veh.mu_x+veh.sens_x*(veh.mu_x_M*g-abs(veh.fz_tyre))).*abs(veh.fz_tyre) ;
        

        Vehicles{r,c} = veh;

    end
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vehicle] = writerows(r,c)
    flag = 0;
    for row = 1:1:size(parameters,1)
        fieldname = fieldnames(row);
        fieldname = strrep(fieldname,' ','');

        if fieldname == 'nf' continue; end %do not write a blank row
        
        %%String or Path fields, just set the value to min and move on
        if ismember(fieldname,stringfields)
                vehicle.(fieldname) = string(parameters(row,2));
            continue %break on every string field
        end

        
        %Define Current min & max values
        min = mincol(row);
        max = maxcol(row);
        if max == 0
            %Set value to min.
            vehicle.(fieldname) = min;
            continue %break on every non sweeping field.
        elseif max < min
            %Throw error if max<min. 
            error('Check input excel sheet \n%s',...
                   fieldname+' max field is less than min field.');
        end
        %increment flag to know which field to iterate.
        flag = flag +1;
        %Set field to a value rising each iteration.
        if flag == 1
            vehicle.(fieldname) = min+((max-min)/(swpits-1))*(r-1);
            vehicle.SweptParameter1 = fieldname;
            if or(fieldname == 'LiftCoefficientCL', fieldname == 'DragCoefficientCD');
                vehicle.SweptParameter1Data = min:((max-min)/(swpits-1)):max;
                vehicle.SweptParameter1Data = flip(vehicle.SweptParameter1Data);
            else
                vehicle.SweptParameter1Data = min:((max-min)/(swpits-1)):max;
            end
            vehicle.SweptParameter2 = "none";
            vehicle.SweptParameter2Data = 0;
        elseif flag == 2
            vehicle.(fieldname) = min+((max-min)/(swpits-1))*(c-1);
            vehicle.SweptParameter2 = fieldname;
            if or(fieldname == 'LiftCoefficientCL',fieldname == 'DragCoefficientCD');
                vehicle.SweptParameter2Data = min:((max-min)/(swpits-1)):max;
                vehicle.SweptParameter2Data = flip(vehicle.SweptParameter2Data);
            else
                vehicle.SweptParameter2Data = min:((max-min)/(swpits-1)):max;
            end
        else
            error('Error in vehicle Array function.')
        end
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = read_torque(workbookFile,sheetName,startRow,endRow)
    % Input handling
    % If no sheet is specified, read first sheet
    if nargin == 1 || isempty(sheetName)
        sheetName = 1;
    end
    % If row start and end points are not specified, define defaults
    if nargin <= 3
        startRow = 2;
        endRow = 1000;
    end
    % Setup the Import Options
    opts = spreadsheetImportOptions("NumVariables", 2);
    % Specify sheet and range
    opts.Sheet = sheetName;
    opts.DataRange = "A" + startRow(1) + ":B" + endRow(1);
    % Specify column names and types
    opts.VariableNames = ["EngineSpeed", "Torque"];
    opts.VariableTypes = ["double", "double"];
    % Setup rules for import
    opts.MissingRule = "omitrow";
    opts = setvaropts(opts, [1, 2], "TreatAsMissing", '');
    % Import the data
    data = readtable(workbookFile, opts, "UseExcel", false);
    for idx = 2:length(startRow)
        opts.DataRange = "A" + startRow(idx) + ":B" + endRow(idx);
        tb = readtable(workbookFile, opts, "UseExcel", false);
        data = [data; tb]; %#ok<AGROW>
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = read_parameters(workbookFile,sheetName,startRow,endRow)
    % Input handling
    % If no sheet is specified, read first sheet
    if nargin == 1 || isempty(sheetName)
        sheetName = 1;
    end
    % If row start and end points are not specified, define defaults
    if nargin <= 3
        startRow = 2;
        endRow = 100;
    end
    % Setup the Import Options
    opts = spreadsheetImportOptions("NumVariables", 3);
    % Specify sheet and range
    opts.Sheet = sheetName;
    opts.DataRange = "B" + startRow(1) + ":D" + endRow(1);
    % Specify column names and types
    opts.VariableNames = ["Variable", "Min", "Max"];
    opts.VariableTypes = ["string", "string", "string"];
    % Setup rules for import
    opts.MissingRule = "fill";
    opts = setvaropts(opts, [1, 2], "TreatAsMissing", '','FillValue','nf'); %no field
    % Import the data
    data = readtable(workbookFile, opts, "UseExcel", false);
    for idx = 2:length(startRow)
        opts.DataRange = "A" + startRow(idx) + ":B" + endRow(idx);
        tb = readtable(workbookFile, opts, "UseExcel", false);
        data = [data; tb]; %#ok<AGROW>
    end

end

end
%Read Vehicle Parameters Function
function [Vehicles,VP] = generatevehicle(filename)

%disp("Vehicle Model Generation Started")
%VehicleParameters Array from importparameters
VP = importparameters(filename);

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
		veh.kw = vp.KerbWeight;
		veh.dm1 = vp.DriverMass1;
		veh.dm2 = vp.DriverMass2;
		veh.M1 = vp.KerbWeight+vp.DriverMass1;
		veh.M2 = vp.KerbWeight+vp.DriverMass2;
		veh.df = vp.FrontMassDistribution/100;
        veh.usmf = vp.FrontUnsprungMass;
        veh.usmfcgh = vp.FrontUnsprungCOG/1000;
        veh.usmr = vp.RearUnsprungMass;
        veh.usmrcgh = vp.RearUnsprungCOG/1000;
        veh.sm1 = veh.M1 - vp.FrontUnsprungMass - vp.RearUnsprungMass; %unsprung mass 1
        veh.sm2 = veh.M2 - vp.FrontUnsprungMass - vp.RearUnsprungMass; %unsprung mass 2
        veh.cgh = vp.CoGHeight; %skidpad and acceleration)
        veh.Ws1 = (veh.M1 -2*veh.usmf - 2*veh.usmr)*9.81; %sprung weight
        veh.Ws2 = (veh.M2 -2*veh.usmf - 2*veh.usmr)*9.81;
        
		% wheelbase
		veh.L = vp.Wheelbase/1000;
        %cg vector
        veh.cg = [-veh.L*(1-veh.df),veh.cgh/1000];
        % track width (used in skidpad only)
        veh.twf = vp.FrontTrackWidth/1000;
        veh.twr = vp.RearTrackWidth/1000;
		% steering rack ratio
		veh.rack = vp.SteeringRatio;
		% aerodynamics
		veh.Cl = vp.LiftCoefficient;
		veh.Cd = vp.DragCoefficient;
		veh.factor_Cl = vp.CLScaleMultiplier;
		veh.factor_Cd = vp.CDScaleMultiplier;
		veh.da = vp.FrontAeroDistribution/100;
		veh.A = vp.FrontalArea;
		veh.rho = vp.AirDensity;
		% brakes
        veh.br_disc_d_f = vp.FrontDiscOuterDiameter/1000;
        veh.br_pad_h_f = vp.FrontPadHeight/1000;
        veh.br_pad_mu_f = vp.FrontPadFrictionCoefficient;
        veh.br_nop_f = vp.FrontCaliperNumberofPistons;
        veh.br_pist_d_f = vp.FrontCaliperPistonDiameter/1000;
        veh.br_mast_d_f = vp.FrontMasterCylinderPistonDiameter/1000;
        
        veh.br_disc_d_r = vp.RearDiscOuterDiameter/1000;
        veh.br_pad_h_r = vp.RearPadHeight/1000;
        veh.br_pad_mu_r = vp.RearPadFrictionCoefficient;
        veh.br_nop_r = vp.RearCaliperNumberofPistons;
        veh.br_pist_d_r = vp.RearCaliperPistonDiameter/1000;
        veh.br_mast_d_r = vp.RearMasterCylinderPistonDiameter/1000;
        
        veh.br_ped_r = vp.PedalRatio;
        veh.br_bb = vp.BrakeBias/100;
		
        % tyres, 
        [~,~,raw] = xlsread(filename);
        idx = find(strcmp(raw,'Tyre Name'));
        veh.tyre = raw{idx + length(raw)};
        veh.tyre_radius = vp.TyreRadius/1000;
        veh.Crr = vp.RollingResistance;
        veh.FTSR = vp.FrontTyreSpringRate;
        veh.RTSR = vp.RearTyreSpringRate;
        
		
		
		% engine
        veh.fuel_type = vp.FuelType;
		veh.factor_power = vp.GlobalPowerFactor;
		veh.n_thermal = vp.ThermalEfficiency/100; %decimal
		veh.fuel_LHV = vp.FuelLowerHeatingValue*1000000; %J/kg
		% drivetrain
		veh.drive = vp.DriveType;
		veh.shift_time = vp.GearShiftTime;
		veh.n_primary = vp.PrimaryGearEfficiency/100;
		veh.n_final = vp.FinalGearEfficiency/100;
		veh.n_gearbox = vp.GearboxEfficiency/100;
		veh.ratio_primary = vp.PrimaryGearReduction;
		veh.ratio_final = vp.FinalGearReduction;
		veh.ratio_gearbox(1) = vp.GearRatio1;
		veh.ratio_gearbox(2) = vp.GearRatio2;
		veh.ratio_gearbox(3) = vp.GearRatio3;
		veh.ratio_gearbox(4) = vp.GearRatio4;
		veh.ratio_gearbox(5) = vp.GearRatio5;
		veh.ratio_gearbox(6) = vp.GearRatio6;
		veh.ratio_gearbox(7) = vp.GearRatio7;
		veh.ratio_gearbox(8) = vp.GearRatio8;
		veh.ratio_gearbox(9) = vp.GearRatio9;
		veh.ratio_gearbox(10) = vp.GearRatio10;
        veh.Rt = vp.RegenTorque;
		for g = 1:length(veh.ratio_gearbox)
			if veh.ratio_gearbox(g) == 0
				veh.nog = g-1;
				veh.ratio_gearbox = veh.ratio_gearbox(veh.ratio_gearbox ~= 0);
				break
			end
		end
		veh.nog = length(veh.ratio_gearbox) ;

% 		%% Brake Model
% 
% 		veh.br_pist_a = veh.br_nop*pi*(veh.br_pist_d/1000)^2/4 ; % [m2]
% 		veh.br_mast_a = pi*(veh.br_mast_d/1000)^2/4 ; % [m2]
% 		veh.beta = veh.tyre_radius/(veh.br_disc_d/2-veh.br_pad_h/2)/veh.br_pist_a/veh.br_pad_mu/4 ; % [Pa/N] per wheel
% 		veh.phi = veh.br_mast_a/veh.br_ped_r*2 ; % [-] for both systems

        %% suspension (used in skidpad calcs only)
        veh.rcf = vp.FrontRollCentreHeight;
        veh.rcr = vp.RearRollCentreHeight;
        veh.kphif = vp.FrontRollStiffness;
        veh.kphir = vp.RearRollStiffness;
        veh.FRHs = 35; %weird fucky error 
        veh.RRHs = 35;
        veh.FSR = vp.FrontSpringRate*175.126835;
        veh.RSR = vp.RearSpringRate*175.126835;
        veh.FMR = vp.FrontMotionRatio;
        veh.RMR = vp.RearMotionRatio;
        veh.FARB = vp.FrontARBRate;
        veh.RARB = vp.RearARBRate;
        veh.hs = perp_dist(veh.rcf/1000,veh.rcr/1000,[veh.L*(1-veh.df),veh.cgh],veh.L);

		%% Steering Model

% 		veh.a = (1-veh.df)*veh.L ; % distance of front axle from center of mass [mm]
% 		veh.b = -veh.df*veh.L ; % distance of rear axle from center of mass [mm]
% 		veh.C = 2*[veh.CF,veh.CF+veh.CR;veh.CF*veh.a,veh.CF*veh.a+veh.CR*veh.b] ; % steering model matrix

		%% Driveline Model
        if veh.fuel_type == "Accumulator" %modified as to have only 1 gear ratio for an EV
            veh.nog = 1;
            veh.ratio_gearbox = veh.ratio_gearbox(1);
        else
        end
        
     
            
		% fetching engine curves
		veh.hpen_speed_curve = vp.hpEngineSpeed ; % [rpm]
		veh.hpen_torque_curve = vp.hpEngineTorque ; % [N*m]
		veh.hpen_power_curve = veh.hpen_torque_curve.*veh.hpen_speed_curve*2*pi/60 ; % [W]
        veh.lpen_speed_curve = vp.lpEngineSpeed ; % [rpm]
		veh.lpen_torque_curve = vp.lpEngineTorque ; % [N*m]
		veh.lpen_power_curve = veh.lpen_torque_curve.*veh.lpen_speed_curve*2*pi/60 ; % [W]
		% memory preallocation
		% wheel speed per gear for every engine speed value
		veh.wheel_speed_gear = zeros(length(veh.hpen_speed_curve),veh.nog) ;
		% vehicle speed per gear for every engine speed value
		veh.vehicle_speed_gear = zeros(length(veh.hpen_speed_curve),veh.nog) ;
		% wheel torque per gear for every engine speed value
		veh.hpwheel_torque_gear = zeros(length(veh.hpen_torque_curve),veh.nog) ;
        veh.lpwheel_torque_gear = zeros(length(veh.lpen_torque_curve),veh.nog) ;
		% calculating values for each gear and engine speed
		for i=1:veh.nog
			veh.wheel_speed_gear(:,i) = veh.hpen_speed_curve/veh.ratio_primary/veh.ratio_gearbox(i)/veh.ratio_final ;
			veh.vehicle_speed_gear(:,i) = veh.wheel_speed_gear(:,i)*2*pi/60*veh.tyre_radius ;
			veh.hpwheel_torque_gear(:,i) = veh.hpen_torque_curve*veh.ratio_primary*veh.ratio_gearbox(i)*veh.ratio_final*veh.n_primary*veh.n_gearbox*veh.n_final ;
            veh.lpwheel_torque_gear(:,i) = veh.lpen_torque_curve*veh.ratio_primary*veh.ratio_gearbox(i)*veh.ratio_final*veh.n_primary*veh.n_gearbox*veh.n_final ;
		end

		% minimum and maximum vehicle speeds
		veh.v_min = min(veh.vehicle_speed_gear,[],'all') ;
		veh.v_max = max(veh.vehicle_speed_gear,[],'all') ;

		% new speed vector for fine meshing
		dv = 0.5/3.6 ;
		veh.vehicle_speed = linspace(veh.v_min,veh.v_max,(veh.v_max-veh.v_min)/dv) ;
		% memory preallocation
		% gear
		veh.gear = zeros(length(veh.vehicle_speed),1) ;
		% engine tractive force
		veh.fxh_engine = zeros(length(veh.vehicle_speed),1) ;
        veh.fxl_engine = zeros(length(veh.vehicle_speed),1) ;
		% engine tractive force per gear
		veh.fxh = zeros(length(veh.vehicle_speed),veh.nog) ;
        veh.fxl = zeros(length(veh.vehicle_speed),veh.nog) ;
		% optimising gear selection and calculating tractive force
       
		for i=1:length(veh.vehicle_speed)
			% tractive force at each speed in each gear
			for j=1:veh.nog
				veh.fxh(i,j) = interp1(veh.vehicle_speed_gear(:,j),veh.hpwheel_torque_gear(:,j)/veh.tyre_radius,veh.vehicle_speed(i),'linear',0) ;
                veh.fxl(i,j) = interp1(veh.vehicle_speed_gear(:,j),veh.lpwheel_torque_gear(:,j)/veh.tyre_radius,veh.vehicle_speed(i),'linear',0) ;
			end
			% getting maximum tractive force and gear
			[veh.fxh_engine(i),veh.gear(i)] = max(veh.fxh(i,:)) ;
            [veh.fxl_engine(i),veh.gear(i)] = max(veh.fxl(i,:)) ;
		end
		% adding values for 0 speed to vectors for interpolation purposes at low speeds
		veh.vehicle_speed = [0;transpose(veh.vehicle_speed)];
		veh.gear = [veh.gear(1);veh.gear] ;
		veh.fxh_engine = [veh.fxh_engine(1);veh.fxh_engine];
        veh.fxl_engine = [veh.fxl_engine(1);veh.fxl_engine] ;
		% final vectors
		% engine speed
		veh.engine_speed = veh.ratio_final*veh.ratio_gearbox(veh.gear)*veh.ratio_primary.*veh.vehicle_speed/veh.tyre_radius*60/2/pi ;
		% wheel torque
		veh.wheel_torqueh = veh.fxh_engine*veh.tyre_radius ;
        veh.wheel_torquel = veh.fxl_engine*veh.tyre_radius ;
		% engine torque
		veh.engine_torqueh = veh.wheel_torqueh/veh.ratio_final./veh.ratio_gearbox(veh.gear)/veh.ratio_primary/veh.n_primary/veh.n_gearbox/veh.n_final ;
        veh.engine_torquel = veh.wheel_torquel/veh.ratio_final./veh.ratio_gearbox(veh.gear)/veh.ratio_primary/veh.n_primary/veh.n_gearbox/veh.n_final ;
		% engine power
		veh.engine_powerh = veh.engine_torqueh.*veh.engine_speed*2*pi/60 ;
        veh.engine_powerl = veh.engine_torquel.*veh.engine_speed*2*pi/60 ;

		%% Shifting Points and Rev Drops
        
        
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
% 		% Z axis
% 		veh.fz_mass = -veh.M1*g ;
% 		veh.fz_aero = 1/2*veh.rho*veh.Cl*veh.A*veh.vehicle_speed.^2 ;
% 		veh.fz_total = veh.fz_mass+veh.fz_aero ;
% 		veh.fz_tyre = (veh.factor_drive*veh.fz_mass+veh.factor_aero*veh.fz_aero)/veh.driven_wheels ;
% 		% x axis
% 		veh.fx_aero = 1/2*veh.rho*veh.Cd*veh.A*veh.vehicle_speed.^2 ;
% 		veh.fx_roll = veh.Cr*abs(veh.fz_total) ;
% 		veh.fx_tyre = veh.driven_wheels*(veh.mu_x+veh.sens_x*(veh.mu_x_M*g-abs(veh.fz_tyre))).*abs(veh.fz_tyre) ;
		
		
      %% VD Calcultions %%
      veh.Fwr = (veh.FTSR*(veh.FSR/veh.FMR^2))/(veh.FTSR+(veh.FSR/veh.FMR^2)); %front wheel rate
      veh.Rwr = (veh.RTSR*(veh.RSR/veh.RMR^2))/(veh.RTSR+(veh.RSR/veh.RMR^2)); %rear wheel rate
      veh.FsRR = (pi*veh.Fwr*veh.twf^2)/360; %front spring Roll rate
      veh.RsRR = (pi*veh.Rwr*veh.twr^2)/360;
      veh.kphif = veh.FsRR + veh.FARB; %front roll rate
      veh.kphir = veh.RsRR + veh.RARB; %rear roll rate
      veh.kphi = veh.kphif + veh.kphir; %total roll rate
      veh.RG1 = ((veh.M1*9.81)*veh.hs)/veh.kphi;
      veh.RG2 = ((veh.M2*9.81)*veh.hs)/veh.kphi;
      
      %decoupled%
      
      
      
        
        
        
        Vehicles{r,c} = veh;
        
       

        
    end
end
end








	









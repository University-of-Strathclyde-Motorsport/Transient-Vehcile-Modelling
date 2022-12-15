%% YMD script %%
%this script plots a Yaw Moment Diagram%

clear
clc

addpath("inputs\Tyre")
addpath("inputs\Vehicle")
addpath("Functions")

%% Import Vehcile Model %%
path = pwd;
sheetpath = strcat(path,"\inputs\vehicle\USM23.xlsx");
Vehicles = generatevehicle(sheetpath);
car = Vehicles{1,1};
clear Vehicles
disp('Vehcile Model Loaded Sucesfully')
disp('====================================================================================')

%% Import Tyre Model %%
load("N:\3 - Vehicle Dynamics & Simulation\3 - Laptime Simulation\Oliver Bailey Diss\Lapsim\Inputs\Tyres\Hoosier 43075 16x7.5-10 R25B, 7 inch rim.mat");
disp('Tyre Model Loaded Successfully')
disp('====================================================================================')


%% Setup Input Matrices %%
V = 60/3.6; %60 kph

delta = -50:1:50; %steering angle range, this is assumed that the ackermann is paralell

beta = -10:1:10; %body slip angle

LatAcc = zeros(length(delta),length(beta));
YawRate = LatAcc;


%% Calculate Vehcile Parameters %%

FZ_fs = (car.sm2*9.81*car.df)/2;
FZ_rs = (car.sm2*9.81*(1-car.df))/2;

%% Calculation Loop %%
%guess a beginning lateral velocity and yaw rate
%iterate until conversion

for i = 1:length(delta)
    for j = 1:length(beta)



        for k = 1:1000
            %calculate velocity vectors%
            %the overall velocity are
            Vx = sqrt((V^2)/(1+tand(beta(j)^2))); %longitudinal velocity
            Vy = tand(beta(j))*Vx; %lateral velocity

            if k == 1
                AYG = 0;
                yawrG = 0;
            else
                %lateral acceleration guess
                AYG = Ay/Vx;
                yawrG = yawr;
            end

            %calculate slip angles
            Arf = (Vy + yawrG*a)/(Vx + yawrG*twf/2) - delta(i);
            Alf = (Vy + yawrG*a)/(Vx - yawrG*twf/2) - delta(i);
            Arr = (Vy - yawrG*b)/(Vx + yawrG*twf/2);
            Alr = (Vy - yawrG*b)/(Vx - yawrG*twf/2);


            %calculate aero loads
            DFf = 0.5*car.A*car.rho*car.Cl*V^2*car.da;
            DFr = 0.5*car.A*car.rho*car.Cl*V^2*(1-car.da);
            Drag = 0.5*car.A*car.rho*car.Cd*V^2;

            %calculate lateral load transfer

            LLTf = AYG*(veh.Ws2/(veh.twf))*(veh.hs*((veh.kphif*57.29 + veh.Ws2*veh.hs*(1-veh.df))/(veh.kphif*57.29+veh.kphir*57.29-veh.Ws2*veh.hs)) + (1-veh.df)*(veh.rcf/1000)) + (veh.usmf*9.81*veh.usmfcgh)/(veh.twf);
            LLTr = AYG*(veh.Ws2/(veh.twr))*(veh.hs*((veh.kphir*57.29 + veh.Ws2*veh.hs*veh.df)/(veh.kphif*57.29+veh.kphir*57.29-veh.Ws2*veh.hs)) + veh.df*(veh.rcr/1000)) + (veh.usmf*9.81*veh.usmrcgh)/(veh.twr);


            %calculate wheel loads
            FZfl = -FZ_fs + DFf/2 + LLTf;  %(N)
            FZfr = -FZ_fs + DFf/2 - LLTf; %(N)
            FZrl = -FZ_rs + DFr/2 + LLTr; %(N)
            FZrr = -FZ_rs + DFr/2 - LLTr; %(N)



            %tyre model
            %fronts only need lateral acceleration can be calculated
            %rears need lateral and longitudinal forces
            %it currentyl assumed that the differential is locked

            [FYfl] = TyreModel(Fzfl,0,"Lateral");
            [FYfr] = TyreModel(Fzfr,0,"Lateral");
            [FYrl,FXrl] = TyreModel(Fzrl,Drag/2,"Combined");
            [FYrr,FXrr] = TyreModel(Fzrr,Drag/2,"Combined");

            %resultant forces and moments
            fyf = FYfl + FYfr;
            fyr = FYrl + FYrr;


            %calculate yaw rate and lateral acceleration

            alpha = (fyf*a - FYr*b)/Iz;
            Ay = (fyf+fyr)/car.sm2;
            yawr = Ay/Vx; %yaw rate


            %add in loop ending conditions
            if (Ay > 0.999*AYG) && (Ay < 1.001*AYG) && (yawrG > 0.999*yawrG) && (yawr < 1.001*yawrG)
                break

            end

            

        end
        %end iteration loop

        %once the loop has iterated and the accelerations
















    end
end


%% Post-Processing and Plotting
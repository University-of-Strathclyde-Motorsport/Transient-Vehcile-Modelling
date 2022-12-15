function [ellipse,forces,alpha,SR] = TyreModel(Tyre,FZ,type)

if type == "Combined"

        slips = -0.16:0.01:0.15;
        alpha = -13:1:13;
        
        
        %% initialise arrays %%
        %for solve speed
        
        
        SRN = zeros(length(slips),1);
        ModelFXm = zeros(length(slips),1);
        FXm = zeros(length(slips),1);
        alphaN = zeros(length(alpha),1);
        ModelFY = zeros(length(alpha),1);
        FYm = zeros(length(alpha),1);
        FY = zeros(length(alpha),length(slips));
        FX = zeros(length(alpha),length(slips));
        
        
        
        
        %% expand the tyre model for each load
       
        zerox = find(slips == 0); %index at which zero SR is found.
        
        nlslips = -slips(1:zerox-1); %negative slip ratios (lsip ratio must be positive in equation)
        plslips = slips(zerox+1:end); %positive slip ratios
        
        %% Slip Ratio Modifications %%
        %the negative slip ratios are defined the same as the SAE, so the positive
        %values must be recaclulated%
        plslips = 1  - 1./(plslips + 1);
        
        %% Re Construct slip ratio vector %%
        MNCslips = [nlslips,0,plslips];
        
        Mux = polyval(Tyre.Mux,FZ(i));
        kx = polyval(Tyre.kx,FZ(i));
        
        for slip = 1:length(slips)
            %unormalise slip angles
            SRN(slip) = (kx*slips(slip))/(FZ(i)*Mux);
            
            
            %reference normalised Tyre model Tyre model
            ModelFXm(slip) = Pacejka4(Tyre.CoeffsFX,SRN(slip));
            
            
            %unormalise Tyre forces
            FXm(slip) = ModelFXm(slip)*Mux*FZ(i);
        end
        FXmf = FXm;
        FXm = abs(FXm);
        
        
        zeroy = find(alpha == 0);
        
        Muy = polyval(Tyre.Muy,FZ(i));
        CS = polyval(Tyre.CS,FZ(i));
        
        
        for sa = 1:length(alpha)
            %unormalise slip angles
            alphaN(sa) = (CS.*tand(alpha(sa)))*57.29/(FZ(i)*Muy);
            
            
            %reference normalised tyre model tyre model
            ModelFY(sa) = Pacejka4(Tyre.CoeffsFY,alphaN(sa));
            
            
            %unormalise tyre forces
            FYm(sa) = ModelFY(sa)*Muy*FZ(i);
            
            
        end
        FYmf = FYm;
        FYm = abs(FYm);
        alpha = abs(alpha);
        CS = CS*57.29;
        
        %% input the nicolas comstock model %%
        %the slip angle is known so all that needs to be done is to iterate
        %the slip ratio value until the FX matches that of the 

                                
                FY(j,k) = (FXm(k) * FYm(j)  / sqrt(MNCslips(k)^2 * FYm(j)^2 + ((FXm(k)^2)*(tand(alpha(j)))^2))) * (sqrt((1-MNCslips(k))^2 * cosd(alpha(j))^2 * FYm(j)^2 + sind(alpha(j))^2 * kx^2) /(kx * cosd(alpha(j))));
                FX(j,k) = (FXm(k) * FYm(j) * MNCslips(k))  / sqrt(MNCslips(k)^2 * FYm(j)^2 + ((FXm(k)^2)*(tand(alpha(j)))^2)) * (sqrt(MNCslips(k)^2 * CS^2 +  (1 - MNCslips(k))^2 * cosd(alpha(j))^2 * FXm(k)^2)/(MNCslips(k)*CS));
                

        
        
   
end
    
    
    
end
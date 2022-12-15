%Read Vehicle Parameters Function
% - Murray
%This function creates an array of base parameters, no calculations or model generation is done here.
%The car 'model' is created in generatevehicle.
function [Vehicles,dispparameters] = importparameters(filename)
%%read the sheet as strings as to be able to import the excel file into the
%app
dispparameters = string(readcell(filename,'Sheet',1,'Range','B2:E62'));



%% Read vehicle file
parameters = table2cell(read_parameters(filename,'Parameters'));
hptorque = read_torque(filename,'High Power Map');
lptorque = read_torque(filename,'Low Power Map');
swpits = str2double(string(parameters(3,2)));

%Create new arrays for each column. Not nice but matlab 
%wasn't playing well when working directly with parameters().
fieldnames = string(parameters(:,1));
mincol = str2double(string(parameters(:,2)));
maxcol = string(parameters(:,3));
maxcol = fillmissing(maxcol,'constant',"0");
maxcol = str2double(maxcol);
stringfields = ["Name","Type","Aeromap","FrontTyreModel","RearTyreModel","FrontKinematicsData","RearKinematicsData","FuelType","DriveType"];

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
        Vehicle.hpEngineSpeed = table2array(hptorque(:,1));
        Vehicle.hpEngineTorque = table2array(hptorque(:,2));
        Vehicle.lpEngineSpeed = table2array(lptorque(:,1));
        Vehicle.lpEngineTorque = table2array(lptorque(:,2));
        Vehicles{r,c} = Vehicle;
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
%Labbook JN 59-61
% the purpose of this file is to assign the measurement performed in your
% lab to the correct "physical" meaning. The important point is, that the
% filename that is built by this function points to the correct file! Our
% file names look like: MonDY_num.tim So a 3 letter code for the month,
% followed by a 2 number code for the day followed by an underscore and
% then a number of the measurement, starting at 1 ending at 999... 
function [temprangeall,monat_messung,tag_messung,name_messung_air,name_messung_reference,name_messung_sample,...
    ordner,measurementtype,ordnerspeichern,nummer_messung_sample,nummer_messung_reference,nummer_messung_air,probenname] ...
    = Assigment_sheet_TDS(nummerderreihe,temp)

global dteflon
global dpila

%% these two values are only needed for TDS_Pellet 
global mischung % mass of the sample material
global volumenfaktor %volume fraction
%% 
measurementtype = 'TDS_results\'; %sample name for folder creation
ordnerspeichern = 'C:\ExperimentalData_backuped\jens\pila2019\'; %output folder
temprangeall =[3301]; % You can define here sets of 
%% this block is needed to ensure a smooth running of the assigment file with TDS_Pellets. Originally each main program hat its own assigment_sheet.m syntax ... 
tag_messung =666;monat_messung ='this is a dummy';
nummer_messung_reference =666;nummer_messung_air = 667; nummer_messung_sample =665;
ordner ='dummy';probenname ='dummy';

if nummerderreihe ==1 % Supports multiple series, so one file can be used for multiple samples    
    % Historically, temp is used as unique identifier for each measurement,
    % see the read.me file 
    % Each measurement consits of 2 or 3 individual measurements! The first    
    %% This block must be defined for each measurement (consiting of Sample/Referance/Air measurement)
    if temp ==3301 % Sample 3, Measurement 1
        nummer_messung_sample = 11; %this is used later to build the correct file name ! 
        nummer_messung_reference = 11; %this is used later to build the correct file name ! 
        nummer_messung_air = 10; %this is used later to build the correct file name !
    end
    if temp ==3302 % Sample 3, Measurement 2
        nummer_messung_sample = 12;
        nummer_messung_reference = 12;
        nummer_messung_air = 10;
    end
    
    if temp ==5301 % Sample 5, Measurement 1
        nummer_messung_sample = 3;
        nummer_messung_reference = 3;
        nummer_messung_air = 2;
    end
    %% This block must be defined once per sample
    if temp > 3000 %sample 3
        monat_messung = 'Jan';
        tag_messung = 30; %day of the measurement
        probenname ='some_name_to_ID_the sample'; %sample name (not to be confused with file name!)
        dteflon =4947.9; % The thickness of the Substrate
        measurementtype = '';
        dpila =10; %the thickness of the sample material
        %input folder
        ordner= 'C:\ExperimentalData_backuped\Jens\pila2019\emptyholder\'; %Here, change me! Change me!
        % this must point to an existing folder!
    end
    if temp > 5000 %sample 5
        monat_messung = 'Jan';
        tag_messung = 30; %day of the measurement
        probenname ='some_name_to_ID_the sample'; %sample name (not to be confused with file name!)
        dteflon =4947.9; % The thickness of the Substrate
        measurementtype = '';
        dpila =10; %the thickness of the sample material
        %input folder
        ordner= 'C:\ExperimentalData_backuped\Jens\pila2019\emptyholder\'; %Here, change me! Change me!
        % this must point to an existing folder!
    end
    %% This block builts the correct file name from the previous input info. You will have to edit this part massively!
    % Check with the file names create. It might be easier to hard-code the
    % names of the measurements instead of using this if/then/else
    % structure. Our measurements are named by this scheme, but your's
    % maybe not ... So you have to decide what is more convenient for you
    if tag_messung >9
        if nummer_messung_reference>9
            name_messung_reference = [ordner,measurementtype,monat_messung,num2str(tag_messung),'_0',num2str(nummer_messung_reference),'.tim'];
        else
            name_messung_reference = [ordner,measurementtype,monat_messung,num2str(tag_messung),'_00',num2str(nummer_messung_reference),'.tim'];
        end
        if nummer_messung_reference>99
            name_messung_reference = [ordner,measurementtype,monat_messung,num2str(tag_messung),'_',num2str(nummer_messung_reference),'.tim'];
        end
    else
        if nummer_messung_reference>9
            name_messung_reference = [ordner,measurementtype,monat_messung,'0',num2str(tag_messung),'_0',num2str(nummer_messung_reference),'.tim'];
        else
            name_messung_reference = [ordner,measurementtype,monat_messung,'0',num2str(tag_messung),'_00',num2str(nummer_messung_reference),'.tim'];
        end
        if nummer_messung_reference>99
            name_messung_reference = [ordner,measurementtype,monat_messung,'0',num2str(tag_messung),'_',num2str(nummer_messung_reference),'.tim'];
        end
    end
    if tag_messung >9
        if nummer_messung_sample>9
            name_messung_sample = [ordner,measurementtype,monat_messung,num2str(tag_messung),'_0',num2str(nummer_messung_sample),'.tim'];
        else
            name_messung_sample = [ordner,measurementtype,monat_messung,num2str(tag_messung),'_00',num2str(nummer_messung_sample),'.tim'];
        end
        if nummer_messung_sample>99
            name_messung_sample = [ordner,measurementtype,monat_messung,num2str(tag_messung),'_',num2str(nummer_messung_sample),'.tim'];
        end
    else
        if nummer_messung_sample>9
            name_messung_sample = [ordner,measurementtype,monat_messung,'0',num2str(tag_messung),'_0',num2str(nummer_messung_sample),'.tim'];
        else
            name_messung_sample= [ordner,measurementtype,monat_messung,'0',num2str(tag_messung),'_00',num2str(nummer_messung_sample),'.tim'];
        end
        if nummer_messung_sample>99
            name_messung_sample = [ordner,measurementtype,monat_messung,'0',num2str(tag_messung),'_',num2str(nummer_messung_sample),'.tim'];
        end
    end
    if tag_messung >9
        if nummer_messung_air>9
            name_messung_air = [ordner,measurementtype,monat_messung,num2str(tag_messung),'_0',num2str(nummer_messung_air),'.tim'];
        else
            name_messung_air = [ordner,measurementtype,monat_messung,num2str(tag_messung),'_00',num2str(nummer_messung_air),'.tim'];
        end
        if nummer_messung_air>99
            name_messung_air = [ordner,measurementtype,monat_messung,num2str(tag_messung),'_',num2str(nummer_messung_air),'.tim'];
        end
    else
        if nummer_messung_air>9
            name_messung_air = [ordner,measurementtype,monat_messung,'0',num2str(tag_messung),'_0',num2str(nummer_messung_air),'.tim'];
        else
            name_messung_air= [ordner,measurementtype,monat_messung,'0',num2str(tag_messung),'_00',num2str(nummer_messung_air),'.tim'];
        end
        if nummer_messung_air>99
            name_messung_air= [ordner,measurementtype,monat_messung,'0',num2str(tag_messung),'_',num2str(nummer_messung_air),'.tim'];
        end
    end
    %%
    %% this part is only neede for TDS_Pellet
    mischung(1) = 20;%
mischung(2) = 20;%
massereference =1000;
rho = 1067;%kg/m3 http://www.guidechem.com/dictionary/en/760-78-1.html
rhoishteflon = 3.55; %3.34mm/g = 3.34 mum/mg, email mike williams, 8/17/2016
% here we need the curve from the calibration and the mass offset!
general_mass_offset = 126.6; %determined from fit on calibrated data
    if temp >1000
        mischung(1) = 1126;% check values
        mischung(2) = 38;% check values
        massereference = 1131; %check values
        probenname = 'sample_name';
        rholoc = rho; % the density of the sample material
        tag_messung = 1;
        monat_messung = 'Mar';
    end 
    
     %% Calibration measurements
    massereference = massereference -general_mass_offset;
    percentage_offset = 1-(general_mass_offset/(mischung(1)+mischung(2)));
    mischung(1) = mischung(1)*percentage_offset;
    mischung(2) = mischung(2)*percentage_offset;
    
    radius = (12.79e-3)/2;% in m
    %% calibrated values
    V2 = (mischung(2)*1e-6)/rho; %volume of amino acid
    dzusatz = (V2/(pi*radius*radius))*1e6; %corresponding thickness of a pure amino acid layer
    dteflon = massereference*rhoishteflon; %thickness of the reference teflon pellet
    dDLN =  mischung(1)*rhoishteflon+dzusatz; %thickness of the Aminoacid mixed pellet, calculated from Teflon and amino acid masses
    rho = rho;%kg/m3 http://www.guidechem.com/dictionary/en/760-78-1.html
    r = radius;%10.1039/C1CP20594C mike williams (phys chem chem phys)
    Vtotal = pi*r*r*dDLN*1e-6;%m3
    VDLN = (mischung(2)*1e-6)/rho;
    volumenfaktor = VDLN/Vtotal; %For EMT
    
%% 
end


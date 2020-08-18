%Labbook JN 59-61          ,145
% All rights reserved (C) Jens Neu 2019, the software can be used according to GPL-3 licence.
% The software is provided as is without any warranty or reliability on my
% site. Please report bugs to jens.neu@yale.edu. Please do not distribute.
clear all
close all
symbole = '*';        LW =2; %linewidth for plots
c0 = 2.99792458e8; %speed of light, do not change please ...
%% pointer to data
ordnerspeichern = 'C:\ExperimentalData_backuped\Jens\DLNorleucine_Doped\JN\'; %folder in which the results are stored/saved/exported
%% program task selection
berserker = true; %mark if program runs on berserker or on laptop
evaluatingone =1; %calculating the spectral data, or plotting the spectral data
evaluatall = true; %for one temp or for all
alletemperaturen =true; %all temp or only selected
%% define how much output you want
zwischenplot =0; %showing the individual results, creates a large variety of plots ...
zwischenplottime =1;% timedomainplot
zwischenplotfitten =1; %especially in single evaluation interesting, plots all calculated results (E/Eref,E/Eair...) for all temp. in individual window
%%

gaussfenster =true; %window funktion, if set to true a Gaussian window is used, if false a rect. window is used
fittend =1; % fitts the refractive index of the substrate with a model (debye, cauchy, or Sellmeier), needs fitten! It will also fit the thickness of the sample!
fittenn =1; % optimization algorithm for refractive index
global datenfensternfuerfit
datenfensternfuerfit =1; %gausswindow to kill 1st Teflon/sample reflex doesn't work yet without
nplot =1; %plotting the refractive index
offset = 0; %substract a linear offset, slope is calculated based on first 40 and last 40 points, not tested!
mitdergroessespielen =0; %test IMG and DMG on calculated data, was never properly implemented!

real_imag_seperation =1; %Boolean, true seperates real and imaginary part calculation, false calculates purely numerically and doesn't work yet ...
averaging_in_time =0; %pick wheter time or frequency domain is used for averagingn
howmanyiterations =5;%needs to be defined if averaging in time is 0
ploterrorbar =1; %only works together with averaging in time =0;

singletemp = 10160; %needed if evaluatall =0, to pick the temp of interest
lorentzfit =0; % implemented but not tested; fits a Lorentz function on the peaks
tplus_onair =0; %if time start for tref and tair is not the same, this will adjust the time axis of the air measurement
global sellerie
sellerie =1;% pick cauchy or Sellmeier, Cauchy is not yet implemented!, needs fittend! %please note, that fittend=true , and  sellerie =false, means that the %thickness is fitted with the Cauchy eq. instead of the sellmeier!
nummerderreihe =1; %pick a series that should be processed
global frequenzrange
global freq_for_d_start
global freq_for_d_end
global frequenzschritt
%% frequency range, for plotting, fitting and stuff ...
freq_for_d_start = 0.15; %not sure for what these are used .....
freq_for_d_end = 4.2;    %not sure for what these are used .....
freq_min_fit = freq_for_d_start; %beginning of the calculation region; frequency
freq_max_fit = freq_for_d_end;   %end of the freq data for calculation
freqmin = freq_for_d_start;      %used for most plotting commands. Should be larger or equal to freq_min_fit
freqmax = freq_for_d_end;        %used for most plotting commands. Should be smaller or equal to freq_max_fit
frequenzschritt = 5;             %stepsize in GHz for calculation n,k, and alpha. If to small, calculation is noisy; if to large small features are gone....
ffmin = freq_for_d_start;        %plotting range
ffmax =3.2;                      %plotting range
%% end input
%% Read out data from zuordnungmessung... and define frequency steps
if datenfensternfuerfit
    weite =8.25; %filter out etalon from teflon window and sample, width for Gaussian window in ps
else
    weite =100; %filter out etalon from antenna
end
if alletemperaturen %this defines whether temprange is defined in the main program (else case) or if the assigment sheet program tells the main program what data it has
    temp=1;
    [temprange,~,~,~,~,~,~,~,~,~,~,~,probenname]= Assigment_sheet_TDS(nummerderreihe,temp); %Please note that this is also used in line 100 and 834, so make sure that you change it here and at line 100 and 834!
else
    nummerderprobe =10;
    temprange1 = [250 200 195 190 185 180 175 170 165 160 155 150 140 130 110]+1000;
    temprange2 = [250 200 195 190 185 180 175 170 165 160 155 150 140 130 110]+10000;
    temprange = [temprange1 temprange2];
    temprange=temprange2;
    temprange=[10110];
end
%% Set fit and calculation range
frequenzrange = freq_min_fit:(frequenzschritt/1000):freq_max_fit;
omega = 2*pi*frequenzrange*1e12;
global frequenzmitteln
frequenzmitteln =1; %averaging data in frequency
legdichhin{1}='leer'; %initilizes an array that will store the legend entries...
global volumenfaktor %for effective medium theory


%% control variable for fit

global teflon
teflon =0; %see read me
%%
farbe = hsv(length(temprange)); %defines a color map for the plot
%% Body
if evaluatingone
    temp = singletemp;
    if evaluatall
    else
        temprange = temp;
    end
    
    for lauflauf =1:length(temprange);
        temp = temprange(lauflauf) %the current temperature
        
        %% assignmend of measured data
        [~,monat_messung,tag_messung,name_messung_air,name_messung_reference,name_messung_sample,...
            ordner,measurementtype,ordnerspeichern,nummer_messung_sample,nummer_messung_reference,nummer_messung_air,probenname]=...
        Assigment_sheet_TDS(nummerderreihe,temp);
        if averaging_in_time
            zustatz ='TD';
        else
            zustatz ='FD';
        end
        if real_imag_seperation
            zz ='\';
        else
            zz='\numerical\';
        end
        %% transfer to fit program
        global temperatur
        temperatur = temp;
        global n_sample
        n_sample = nummer_messung_sample;
        global n_teflon
        n_teflon = nummer_messung_reference;
        global n_air
        n_air = nummer_messung_air;
        global ordnername
        ordnername = [ordner,measurementtype,monat_messung];
        global tagmessung
        tagmessung = tag_messung;
        global dDLN
        global dteflon
        if n_teflon ==666; %Use this value in the assigment sheet if only air, not air and reference was measured.
            name_messung_reference = name_messung_air;
            normalizedon = 'air';
        else
            normalizedon = 'teflon';
        end
        
        %% assigning data
        clear t_air
        clear A_air
        clear t_reference
        clear A_reference
        clear t_sample
        clear A_sample
        if averaging_in_time
            [t_air, A_air] = Import_for_timfiles(name_messung_air);
            [t_reference, A_reference] = Import_for_timfiles(name_messung_reference);
            [t_sample, A_sample] = Import_for_timfiles(name_messung_sample);
            t_air = t_air +tplus_onair;
        else
            Air_complete = Import_for_timsfiles2([name_messung_air,'s'],howmanyiterations);
            t_air = Air_complete(:,1)+tplus_onair;
            A_air = Air_complete(:,2:end);
            ref_complete = Import_for_timsfiles2([name_messung_reference,'s'],howmanyiterations);
            t_reference = ref_complete(:,1);
            A_reference = ref_complete(:,2:end);
            sample_complete = Import_for_timsfiles2([name_messung_sample,'s'],howmanyiterations);
            t_sample = sample_complete(:,1);
            A_sample = sample_complete(:,2:end);
        end
        tchange = 90;
        if offset
            minimum =mean(mean(A_reference(2:40,:)));
            maximum = mean(mean(A_reference(end-40:end,:)));
            tmin = mean(t_reference(2:40));
            tmax = mean(t_reference(41:80,:));
            m = (maximum-minimum)/(tmax - tmin);
            b = maximum-(m*tmax);
            if averaging_in_time
                Akorrektur_air = t_air*m+b;
            else
                Akorrektur_air = t_air*m+b;
                for laufnochmal=1:min(size(A_air))
                    A_air(:,laufnochmal) = A_air(:,laufnochmal)-Akorrektur_air;
                    A_sample(:,laufnochmal) = A_sample(:,laufnochmal)-Akorrektur_air;
                    A_reference(:,laufnochmal) = A_reference(:,laufnochmal)-Akorrektur_air;
                end
            end
            [m b]
        else
        end
        for lauf=1:length(t_sample)
            if t_sample(lauf) > tchange
                sprungfenster_sample(lauf,1) =0;
            else
                sprungfenster_sample(lauf,1) =1;
            end
        end
        tchange = 90;
        for lauf=1:length(t_air)
            if t_air(lauf) > tchange
                sprungfenster_air(lauf,1) =0;
            else
                sprungfenster_air(lauf,1) =1;
            end
        end
        global freq_air
        global freq_sample
        global freq_reference
        [maxi index] = max(A_sample(:,1));
        deltat =t_sample(index);
        clear fenster
        if gaussfenster
            fenster = exp(-((t_sample-deltat)/weite).^2);
        else
            fenster = sprungfenster_sample;
        end
        clear A_reference_windowed
        clear A_sample_windowed
        numberofiterations = max(size(A_reference(1,:)));
        for argh =1:numberofiterations
            A_reference_windowed(:,argh) = A_reference(:,argh).*fenster;
            A_sample_windowed(:,argh) = A_sample(:,argh).*fenster;
        end
        %% some fft voodoo
        NNT = 2^18;% better 15 or 16, Zero padding!
        freq_reference = faxis(t_reference,NNT);
        freq_sample = faxis(t_sample,NNT);
        freq_air = faxis(t_air,NNT);
        if averaging_in_time
            spec_reference = abs(fft(A_reference,NNT));
            fft_reference = fft(A_reference,NNT);
            spec_reference_windowed = abs(fft(A_reference_windowed,NNT));
            fft_reference_windowed = fft(A_reference_windowed,NNT);
            spec_sample = abs(fft(A_sample,NNT));
            spec_sample_windowed = abs(fft(A_sample_windowed,NNT));
            fft_sample = fft(A_sample,NNT);
            fft_sample_windowed = fft(A_sample_windowed,NNT);
            fft_air = fft(A_air,NNT);
        else
            spec_referencek = abs(fft(A_reference,NNT));
            spec_reference = mean(spec_referencek,2);
            fft_referencek = fft(A_reference,NNT);
            fft_reference = mean(fft_referencek,2);
            spec_reference_windowedk = abs(fft(A_reference_windowed,NNT));
            spec_reference_windowed = mean(spec_reference_windowedk,2);
            fft_reference_windowedk = fft(A_reference_windowed,NNT);
            fft_reference_windowed = mean(fft_reference_windowedk,2);
            spec_samplek = abs(fft(A_sample,NNT));
            spec_sample = mean(spec_samplek,2);
            spec_sample_windowedk = abs(fft(A_sample_windowed,NNT));
            spec_sample_windowed = mean(spec_sample_windowedk,2);
            fft_samplek = fft(A_sample,NNT);
            fft_sample = mean(fft_samplek,2);
            fft_sample_windowedk = fft(A_sample_windowed,NNT);
            fft_sample_windowed = mean(fft_sample_windowedk,2);
            fft_airk = fft(A_air,NNT);
            fft_air = mean(fft_airk,2);
            spec_error = std(spec_sample_windowedk,0,2);
            relative_error = (spec_error./spec_sample_windowed)/2;
        end
        [maxi index] = max(A_air(:,1));
        deltat =t_sample(index);
        clear fenster_air
        if gaussfenster
            fenster_air = exp(-((t_air-deltat)/weite).^2);
        else
            fenster_air = sprungfenster_air;
        end
        clear A_air_windowed
        for argh =1:numberofiterations
            A_air_windowed(:,argh) = A_air(:,argh).*fenster_air;
        end
        
        if averaging_in_time
            fft_air_windowed = fft(A_air_windowed,NNT);
        else
            fft_air_windowedk = fft(A_air_windowed,NNT);
            fft_air_windowed = mean(fft_air_windowedk,2);
            testtest = abs(fft_sample_windowedk)./abs(fft_air_windowedk);
            referenced_spec_error = std(testtest,0,2);
            testtestt = abs(fft_sample_windowedk)./abs(fft_reference_windowedk);
            referenced_spec_errorteflon = std(testtestt,0,2);
        end
        spec_reference_on_air = abs(fft_reference)./abs(fft_air);
        spec_reference_on_air_windowed = abs(fft_reference_windowed)./abs(fft_air_windowed);
        spec_sample_on_air = abs(fft_sample)./abs(fft_air);
        spec_sample_on_air_windowed = abs(fft_sample_windowed)./abs(fft_air_windowed);
        %% Preping data for fit function
        global Eair
        global Ereference
        global Esample
        if datenfensternfuerfit
            Esample = fft_sample_windowed;
            Eair = fft_air_windowed;
            Ereference = fft_reference_windowed;
        else
            Esample = fft_sample;
            Eair = fft_air;
            Ereference = fft_reference;
        end
        global S21teflon
        S21teflon = Ereference./Eair;
        global S21DLN
        S21DLN = Esample./Eair;
        global S21teflonphase
        global S21DLNphase
        global S21DLN_on_Teflon_phase
        global S21DLN_on_Teflon
        S21DLN_on_Teflon = Esample./Ereference;
        S21DLN_on_Teflon_phase = unwrap(angle(S21DLN_on_Teflon));
        S21teflonphase = unwrap(angle(S21teflon));
        S21DLNphase = unwrap(angle(S21DLN));
        global mischung
        %% call fit program
        %% fitting the thickness with cauchy or sellmeier model, both is implemented but not tested. Usually I used the mechanical thickness calculated from a calibration curve
        if fittend %fitting thickness of teflon, using a common model for teflon, usually that function is off
            teflon =1;
            if sellerie
                x= [3400 1.25 0.1 0.2 0.001 0.001]; %teflon_thickness, B1, C1, B2, C2, absortion
                lb = [2500 0 0 0 0 0]; ub = [4500 3 1 1 0.5 1];
                options = optimoptions('fmincon');
                options = optimoptions(options,'Display', 'off');
                options = optimoptions(options,'FunValCheck', 'off');
                [x,fval,exitflag,output,lambda,grad,hessian] = ...
                    fmincon(@retrieval_of_d_mark3,x,[],[],[],[],lb,ub,[],options);
                if teflon
                    dteflon =x(1);
                    dteflon_all(lauflauf) = dteflon;
                    f = abs(freq_for_d_start + freq_for_d_end)/2;
                    wellenl = (c0./f)*1e-6;
                    ll = wellenl.*wellenl;
                    llll = ll.*ll;
                    n2 = sqrt(1+ x(2)*ll./(ll-x(3)) + x(4)*ll./(ll-x(5))) + x(6)*1i;
                    nteflon_all(lauflauf) = n2;
                else
                    dDLN =x(1);
                    dDLN_all(lauflauf) =dDLN;
                    f = abs(freq_for_d_start + freq_for_d_end)/2;
                    wellenl = (c0./f)*1e-6;
                    ll = wellenl.*wellenl;
                    llll = ll.*ll;
                    n2 = sqrt(1+ x(2)*ll./(ll-x(3)) + x(4)*ll./(ll-x(5))) + x(6)*1i;
                    nDLN_all(lauflauf) = n2;
                end
                
            else
                x= [3400 1.4 0.01 0.5]; %teflon_thickness, refrective index no dispersion, absorbtion, daepfung, Cauchy1 cauchy2
                lb = [2500 1.3 0 0]; ub = [4500 1.45 1 3];
                options = optimoptions('fmincon');
                options = optimoptions(options,'Display', 'off');
                options = optimoptions(options,'FunValCheck', 'off');
                [x,fval,exitflag,output,lambda,grad,hessian] = ...
                    fmincon(@retrieval_of_d_mark3,x,[],[],[],[],lb,ub,[],options);
                x(5) =0;
                
                if teflon
                    dteflon =x(1);
                    dteflon_all(lauflauf) = dteflon;
                    f = abs(freq_for_d_start + freq_for_d_end)/2;
                    wellenl = (c0./f)*1e-6;
                    ll = wellenl.*wellenl;
                    llll = ll.*ll;
                    n2 = x(2) + x(4)./ll+x(5)./llll;
                    nteflon_all(lauflauf) = n2 + 1i*x(3);
                else
                    dDLN =x(1);
                    dDLN_all(lauflauf) =dDLN;
                    f = abs(freq_for_d_start + freq_for_d_end)/2;
                    wellenl = (c0./f)*1e-6;
                    ll = wellenl.*wellenl;
                    llll = ll.*ll;
                    n2 = x(2) + x(4)./ll+x(5)./llll;
                    nDLN_all(lauflauf) = n2 + 1i*x(3);
                end
            end
            teflon =0;
            if sellerie  %assuming sellmeier eq for the sample pellet, allows to fit thickness and n simulatniously
                x= [3400 1.25 0.0 0.2 0.001 0.001]; %teflon_thickness, B1, C1, B2, C2, absortion
                lb = [2500 0 0 0 0 0]; ub = [4500 3 1 1 0.5 1];
                options = optimoptions('fmincon');
                options = optimoptions(options,'Display', 'off');
                options = optimoptions(options,'FunValCheck', 'off');
                [x,fval,exitflag,output,lambda,grad,hessian] = ...
                    fmincon(@retrieval_of_d_mark3,x,[],[],[],[],lb,ub,[],options);
                if teflon
                    dteflon =x(1);
                    dteflon_all(lauflauf) = dteflon;
                    f = abs(freq_for_d_start + freq_for_d_end)/2;
                    wellenl = (c0./f)*1e-6;
                    ll = wellenl.*wellenl;
                    llll = ll.*ll;
                    n2 = sqrt(1+ x(2)*ll./(ll-x(3)) + x(4)*ll./(ll-x(5))) + x(6)*1i;
                    nteflon_all(lauflauf) = n2;
                else
                    dDLN =x(1);
                    dDLN_all(lauflauf) =dDLN;
                    f = abs(freq_for_d_start + freq_for_d_end)/2;
                    wellenl = (c0./f)*1e-6;
                    ll = wellenl.*wellenl;
                    llll = ll.*ll;
                    n2 = sqrt(1+ x(2)*ll./(ll-x(3)) + x(4)*ll./(ll-x(5))) + x(6)*1i;
                    nDLN_all(lauflauf) = n2;
                end
            else %Thickness is fitted but how ....
                x= [3400 1.3 0.01 0.5]; %teflon_thickness, refrective index no dispersion, absorbtion, daepfung, Cauchy1 cauchy2
                lb = [2500 1.2 0 0]; ub = [4500 1.5 1 3];
                options = optimoptions('fmincon');
                options = optimoptions(options,'Display', 'off');
                options = optimoptions(options,'FunValCheck', 'off');
                [x,fval,exitflag,output,lambda,grad,hessian] = ...
                    fmincon(@retrieval_of_d_mark3,x,[],[],[],[],lb,ub,[],options);
                x(5) =0;
                if teflon
                    dteflon =x(1);
                    dteflon_all(lauflauf) = dteflon;
                    feinzelp = abs(freq_for_d_start + freq_for_d_end)/2;
                    wellenl = (c0./feinzelp)*1e-6;
                    ll = wellenl.*wellenl;
                    llll = ll.*ll;
                    n2 = x(2) + x(4)./ll+x(5)./llll;
                    nteflon_all(lauflauf) = n2 + 1i*x(3);
                else
                    dDLN =x(1);
                    dDLN_all(lauflauf) =dDLN;
                    feinzelp = abs(freq_for_d_start + freq_for_d_end)/2;
                    wellenl = (c0./feinzelp)*1e-6;
                    ll = wellenl.*wellenl;
                    llll = ll.*ll;
                    n2 = x(2) + x(4)./ll+x(5)./llll;
                    nDLN_all(lauflauf) = n2 + 1i*x(3);
                end
            end
            %%plot fitted data
            if zwischenplot
                freqstart =freq_for_d_start;
                freqende = freq_for_d_end;
                freqstartindex = find(freq_air>freqstart,1,'first');
                freqendeindex = find(freq_air>freqende,1,'first');
                f = freq_air(freqstartindex:freqendeindex);
            end
            lege{lauflauf}= [num2str(temprange(lauflauf)),' K '];
            save([ordnerspeichern,num2str(temp),'thicknessteflon'],'dteflon','-ascii');
            save([ordnerspeichern,num2str(temp),'thicknesssample'],'dDLN','-ascii');
        end
        %% fitting the refractive index.
        if fittenn
            global specnorm
            specnorm = max(abs(Eair));
            global phasenorm
            phasenorm = max(unwrap(angle((Eair))));
            global frequenzpunkt
            % seed for n and k
            nteflon(1) = 1.5;
            nDLN(1)  = 1.5;
            nDLN_effectivmedium(1) =1.5;
            nDLN_referenceT(1) =1.5;
            kteflon(1) = 0;
            kDLN(1) =0;
            kDLN_effectivmedium(1) =0;
            kDLN_referenceT(1) =0;
            
            %% paramter for Maxwell-Garnett
            
            for fitschleife =2:length(frequenzrange)
                frequenzpunkt = frequenzrange(fitschleife);
                % refractive index of mixture and reference
                teflon =1;
                [nteflon(fitschleife),kteflon(fitschleife)] = retrieval_of_n_airref_nooptimization(nteflon(fitschleife-1),kteflon(fitschleife-1));
                teflon =0;
                [nDLN(fitschleife),kDLN(fitschleife)] = retrieval_of_n_airref_nooptimization(nDLN(fitschleife-1),kDLN(fitschleife-1));
                %Maxwell Garnett and Bruggeman
                [nDLN_effectivmedium(fitschleife),kDLN_effectivmedium(fitschleife),nDLN_referenceT(fitschleife),kDLN_referenceT(fitschleife),nDLN_effectivemediumT(fitschleife),kDLN_effectivemediumT(fitschleife),...
                    nDLN_effectivmediumB(fitschleife),kDLN_effectivmediumB(fitschleife),nDLN_effectivemediumBT(fitschleife),kDLN_effectivemediumBT(fitschleife)] =...
                    retrieval_of_n_relativtoteflon_nooptimization(nteflon(fitschleife),kteflon(fitschleife),nDLN(fitschleife),kDLN(fitschleife),nDLN_referenceT(fitschleife-1),kDLN_referenceT(fitschleife-1));
                % first pair is MG
            end
            if n_teflon ==666
                kDLN_effectivemediumT = kDLN_effectivmedium;
                kDLN_effectivemediumBT =kDLN_effectivmediumB;
                nDLN_effectivemediumT = nDLN_effectivemedium;
                nDLN_effectivemediumBT = nDLN_effectivemediumB;
            end
            alpha_effektivmedium = (2*kDLN_effectivmedium.*omega/c0)*1e-2; % Maxwell Garnett
            alpha_effektivmediumT = (2*kDLN_effectivemediumT.*omega/c0)*1e-2;%per cm WTF? That is per definition the same as the previous line.... why did I do that this way? Did I forgett to implment the normalizing on Teflon? Or did I kick it out... we will never know...
            alpha_effektivmediumB = (2*kDLN_effectivmediumB.*omega/c0)*1e-2;
            alpha_effektivmediumBT = (2*kDLN_effectivemediumBT.*omega/c0)*1e-2;%per cm
            
            if lorentzfit
                global omegarange
                global epsilon_i
                global epsilon_r
                lorentzrange = [0.5 2];
                indexstart = find(frequenzrange>lorentzrange(1),1,'first');
                indexend = find(frequenzrange>lorentzrange(2),1,'first');
                kDLN_effectivemedium_mean = (kDLN_effectivmedium+kDLN_effectivemediumT+kDLN_effectivmediumB+kDLN_effectivemediumBT)/4;
                
                nDLN_effectivemedium_mean = (nDLN_effectivmedium+nDLN_effectivmediumB+nDLN_effectivemediumBT +nDLN_effectivemediumT)/4;
                epsilon_r = nDLN_effectivemedium_mean(indexstart:indexend).^2- kDLN_effectivemedium_mean(indexstart:indexend).^2;
                epsilon_i = 2* nDLN_effectivemedium_mean(indexstart:indexend).*kDLN_effectivemedium_mean(indexstart:indexend).^2;
                
                omegarange = frequenzrange(indexstart:indexend)*2*pi;
                x = [2.7 2.2 1.5 1];
                lb = [1 0.75 1.3 0]; ub = [5 5 1.6 Inf];
                options = optimoptions('fmincon');
                options = optimoptions(options,'Display', 'off');
                options = optimoptions(options,'FunValCheck', 'off');
                [x,fval,exitflag,output,lambda,grad,hessian] = ...
                    fmincon(@lorentzianfit,x,[],[],[],[],lb,ub,[],options);
                offset = x(1);
                omegap = x(2);
                omega0 = x(3);
                gammar = x(4);
                nenner = (omega0^2-omegarange.^2).^2 + omegarange.^2 *gammar^2;
                zaehlerreal = (omega0^2-omegarange.^2)*omegap^2;
                zaehlerimag = omegarange.*omegap.*omegap.*gammar;
                epsilon_r_l = offset+ (zaehlerreal./nenner);
                epsilon_i_l = zaehlerimag./nenner;
                figure(808+lauflauf)
                hold on
                plot(frequenzrange(indexstart:indexend),epsilon_r,'d','color','b');
                plot(frequenzrange(indexstart:indexend),epsilon_r_l,'color','c');
                plot(frequenzrange(indexstart:indexend),epsilon_i,'+','color','r');
                plot(frequenzrange(indexstart:indexend),epsilon_i_l,':','color','m');
            end
            %% saving kind of everything, don't change is needed for the evaluationon=0 setting
            save([ordnerspeichern,num2str(temp),'refractiveindexsample'],'nDLN','-ascii');
            save([ordnerspeichern,num2str(temp),'refractiveindexteflon'],'nteflon','-ascii');
            save([ordnerspeichern,num2str(temp),'absorptionindexsample'],'kDLN','-ascii');
            save([ordnerspeichern,num2str(temp),'absorptionindexteflon'],'kteflon','-ascii');
            save([ordnerspeichern,num2str(temp),'MG_onT'],'nDLN_effectivemediumT','-ascii');
            save([ordnerspeichern,num2str(temp),'MG_onT_absorb'],'kDLN_effectivemediumT','-ascii');
            save([ordnerspeichern,num2str(temp),'alphaMG_onT'],'alpha_effektivmediumT','-ascii');
            save([ordnerspeichern,num2str(temp),'alphaMG_onT_absorb'],'alpha_effektivmedium','-ascii');
            save([ordnerspeichern,num2str(temp),'B_onair'],'nDLN_effectivmediumB','-ascii');
            save([ordnerspeichern,num2str(temp),'B_onair_absorb'],'kDLN_effectivmediumB','-ascii');
            save([ordnerspeichern,num2str(temp),'B_onT'],'nDLN_effectivemediumBT','-ascii');
            save([ordnerspeichern,num2str(temp),'B_onT_absorb'],'kDLN_effectivemediumBT','-ascii');
            save([ordnerspeichern,num2str(temp),'alphaB_onA'],'alpha_effektivmediumBT','-ascii');
            save([ordnerspeichern,num2str(temp),'alphaB_onT_absorb'],'alpha_effektivmediumB','-ascii');
            save([ordnerspeichern,num2str(temp),'dDLN'],'dDLN','-ascii');
            save([ordnerspeichern,num2str(temp),'dteflon'],'dteflon','-ascii');
            save([ordnerspeichern,num2str(temp),'normalizedon'],'normalizedon','-ascii');
            save([ordnerspeichern,num2str(temp),'MG_onair'],'nDLN_effectivmedium','-ascii');
            save([ordnerspeichern,num2str(temp),'MG_onair_absorb'],'kDLN_effectivmedium','-ascii');
            save([ordnerspeichern,num2str(temp),'freqfit'],'frequenzrange','-ascii')
            save([ordnerspeichern,num2str(temp),'alpha_effektivmedium'],'alpha_effektivmedium','-ascii');
            %% save again the particular interesting results
            save([ordnerspeichern,'pick\',num2str(temp),'_k_teflon'],'kteflon','-ascii');
            save([ordnerspeichern,'pick\',num2str(temp),'_n_teflon'],'nteflon','-ascii');
            save([ordnerspeichern,'pick\',num2str(temp),'_n_peptide'],'nDLN_effectivmedium','-ascii');
            save([ordnerspeichern,'pick\',num2str(temp),'_k_peptide'],'kDLN_effectivmedium','-ascii');
            save([ordnerspeichern,'pick\',num2str(temp),'_frequenz'],'frequenzrange','-ascii')
            save([ordnerspeichern,'pick\',num2str(temp),'_alpha'],'alpha_effektivmedium','-ascii');
            if averaging_in_time
            else %this calculates error bars
                loopvar =1;
                alpha_err(1) =0;
                for fitschleife =2:length(frequenzrange)
                    frequenzpunkt = frequenzrange(fitschleife);
                    loopvar = loopvar+1;
                    index_low = find(frequenzpunkt>(freq_air+frequenzschritt/2000),1,'last');
                    index_high = find(frequenzpunkt>(freq_air-frequenzschritt/2000),1,'last');
                    alpha_err(loopvar) = mean(relative_error(index_low:index_high));
                end
                alpha_err = alpha_err.*alpha_effektivmedium;
                save([ordnerspeichern,'pick\',num2str(temp),'_alpha_error'],'alpha_effektivmedium','-ascii');
            end
            
            %% Lots of Plotting ...
            
            if zwischenplot
                figure(444333)
                hold on
                lege{lauflauf}= [num2str(temp),' K '];
                plot(frequenzrange,nteflon,':','color',farbe(lauflauf,:))
                plot(frequenzrange,nDLN,'color',farbe(lauflauf,:))
                legedichhin{(lauflauf*4)-3}=['ref on air MG ',num2str(temp),' K '];
                legedichhin{(lauflauf*4-2)}=['ref on ',normalizedon,'MG',num2str(temp),' K '];
                legedichhin{(lauflauf*4-1)}=['ref on air B',num2str(temp),' K '];
                legedichhin{(lauflauf*4)}=['ref on ',normalizedon,' B',num2str(temp),' K '];
                xlabel('frequency')
                ylabel('refractive index pellet')
                figure(666)
                hold on
                plot(frequenzrange,nteflon,'color',farbe(lauflauf,:))
                xlabel('frequency [THz]')
                ylabel('refractive index teflon')
                title('teflon')
                
            else
                if nplot
                    if zwischenplotfitten
                        figure(444333)
                        hold on
                        plot(frequenzrange,nteflon,'+','color',farbe(lauflauf,:))
                        figure(444444)
                        hold on
                        plot(frequenzrange,nDLN,'x','color',farbe(lauflauf,:))
                        plot(frequenzrange,nDLN_referenceT,'d','color',farbe(lauflauf,:))
                        figure(444334)
                        hold on
                        plot(frequenzrange,kteflon,'+','color',farbe(lauflauf,:))
                        figure(444443)
                        hold on
                        plot(frequenzrange,kDLN,'x','color',farbe(lauflauf,:))
                        plot(frequenzrange,kDLN_referenceT,'d','color',farbe(lauflauf,:))
                        figure(76)
                        plot(frequenzrange,kteflon,'-','color',farbe(lauflauf,:))
                        hold on
                        plot(frequenzrange,kDLN,'x','color',farbe(lauflauf,:))
                        plot(frequenzrange,kDLN_referenceT,'d','color',farbe(lauflauf,:))
                        title('teflon and sample in the same graph')
                    end
                    figure(4444433)
                    hold on
                    plot(frequenzrange,kDLN_effectivmedium,'x','color',farbe(lauflauf,:))
                    plot(frequenzrange,kDLN_effectivemediumT,'d','color',farbe(lauflauf,:))
                    plot(frequenzrange,kDLN_effectivmediumB,'+','color',farbe(lauflauf,:))
                    plot(frequenzrange,kDLN_effectivemediumBT,'p','color',farbe(lauflauf,:))
                    figure(4444444)
                    hold on
                    plot(frequenzrange,nDLN_effectivmedium,'x','color',farbe(lauflauf,:))
                    plot(frequenzrange,nDLN_effectivemediumT,'d','color',farbe(lauflauf,:))
                    plot(frequenzrange,nDLN_effectivmediumB,'+','color',farbe(lauflauf,:))
                    plot(frequenzrange,nDLN_effectivemediumBT,'p','color',farbe(lauflauf,:))
                    figure(4343)
                    if temp == 10140
                        global alpha_tester
                        alpha_tester = alpha_effektivmedium;
                    end
                    drei=4
                    hold on
                    
                    if averaging_in_time
                        plot(frequenzrange,alpha_effektivmedium,'color',farbe(lauflauf,:))
                        %plot(frequenzrange,alpha_effektivmediumT,'d','color',farbe(lauflauf,:))
                        %plot(frequenzrange,alpha_effektivmediumB,'+','color',farbe(lauflauf,:))
                        %plot(frequenzrange,alpha_effektivmediumBT,'p','color',farbe(lauflauf,:))
                    else
                        errorbar(frequenzrange(1:10:end),alpha_effektivmedium(1:10:end),alpha_err(1:10:end),'color',farbe(lauflauf,:))
                        figure(temp*10+9)
                        errorbar(frequenzrange,alpha_effektivmedium,alpha_err,'color',farbe(lauflauf,:))
                        xlim([ffmin ffmax])
                        ylabel('absorptioncoef. \alpha [cm^{-1}]')
                        xlabel('frequency [THz]')
                        title(['retrieved absorptioncoefficient,',probenname,' EMT, MG'])
                    end
                    
                end
                
                
                tempforplot= temprange(lauflauf);
                zusat ='s0';
                
                zusat = probenname;
                drei =4;
                lege{lauflauf}= [zusat,num2str(tempforplot),' K '];
                legedichhin{(lauflauf*4)-3}=[zusat,'ref on air MG ',num2str(tempforplot),' K '];
                legedichhin{(lauflauf*4-2)}=[zusat,'ref on teflon MG',num2str(tempforplot),' K '];
                legedichhin{(lauflauf*4-1)}=[zusat,'ref on air B',num2str(tempforplot),' K '];
                legedichhin{(lauflauf*4)}=[zusat,'ref on teflon B',num2str(tempforplot),' K '];
                legende2{lauflauf} = [probenname,' at ',num2str(tempforplot),' K '];
                
            end
        end
        %ffmin = 0.3;
        % ffmax =2.6;
        if nplot
            if zwischenplotfitten
                figure(444333)
                box on
                legend(lege);
                xlim([ffmin ffmax])
                title(['retrieved realpart of n, teflon, for thickness',num2str(dteflon)])
                figure(444444)
                box on
                legend(lege);
                xlim([ffmin ffmax])
                title(['retrieved realpart of n,',probenname,' for thickness',num2str(dDLN)])
                xlabel('frequency [THz]')
                figure(444334)
                box on
                legend(lege);
                xlim([ffmin ffmax])
                title(['retrieved imaginarypart of n, teflon, for thickness',num2str(dteflon)])
                xlabel('frequency [THz]')
                figure(444443)
                box on
                legend(lege);
                xlim([ffmin ffmax])
                title(['retrieved imaginarypart of n,',probenname,', for thickness',num2str(dDLN)])
                xlabel('frequency [THz]')
                figure(76)
                box on
                legend(lege);
                xlim([ffmin ffmax])
                title(['Teflon and retrieved imaginarypart of n,',probenname,', for thickness',num2str(dDLN)])
                xlabel('frequency [THz]')
            end
            figure(4444444)
            box on
            legend(legedichhin);
            xlim([ffmin ffmax])
            title(['retrieved realpart of n',probenname,', for thickness ',num2str(dDLN),' effective medium'])
            ylabel('phasecoefficient n^, ')
            xlabel('frequency [THz]')
            figure(4444433)
            box on
            legend(legedichhin);
            xlim([ffmin ffmax])
            title(['retrieved imagpart of n, ',probenname,', for thickness ',num2str(dDLN),' effective medium'])
            ylabel('extinctioncoef. n" ')
            xlabel('frequency [THz]')
            figure(4343)
            box on
            %legend(legende2,'Location','NorthWest');
            xlim([ffmin ffmax])
            ylabel('absorptioncoef. \alpha [cm^{-1}]')
            xlabel('frequency [THz]')
            title(['retrieved absorptioncoefficient,',probenname,' EMT, MG'])
            legend(num2str((mod(temprange,1000)')))
        end
        %% plotting
        if zwischenplottime
            figure(temp*10+1)
            hold on
            plot(t_air,A_air*1000,'k')
            plot(t_reference,A_reference*1000,'b')
            %plot(t_reference,A_reference_windowed,'b')
            plot(t_sample,A_sample*1000,'r')
            plot(t_reference,fenster*1e-2*1000,'g')
            %plot(t_sample,A_sample_windowed,'m')
            % plot(t_sample,fenster*1e-4,'g')
            figure(temp*10+2)
            hold on
            plot(freq_air,abs(fft_air),'c')
            plot(freq_reference,spec_reference,'k')
            plot(freq_reference,spec_reference_windowed,'b')
            xlim([freqmin freqmax])
            plot(freq_sample,spec_sample,'r')
            if averaging_in_time
                plot(freq_sample,spec_sample_windowed,'m')
            else
                index = find(freq_sample>ffmax+1,1);
                h = errorbar(freq_sample(1:10:index),spec_sample_windowed(1:10:index),spec_error(1:10:index),'m')
                h.CapSize = 12;
            end
            
        end
        
        spectral_intensity_absorbtion = spec_sample./spec_reference;
        spectral_intensity_absorbtion_windowed = spec_sample_windowed./spec_reference_windowed;
        
        if zwischenplot
            figure(temp*10+4)
            plot(freq_sample,spectral_intensity_absorbtion,'k')
            hold on
            plot(freq_sample,spectral_intensity_absorbtion_windowed,'b')
            if averaging_in_time
                plot(freq_sample,spec_sample_on_air,'r')
            else
                index = find(freq_sample>ffmax+1,1);
                h = errorbar(freq_sample(1:10:index),spec_sample_on_air(1:10:index),referenced_spec_error(1:10:index),'r')
                h.CapSize = 12;
            end
            xlim([freqmin freqmax])
            ylim([0.1 1.1])
            legend('Sample spectrum divided by teflon','Sample spectrum divided by teflon widowed')
        end
        fft_relative = abs(fft_sample)./abs(fft_reference);
        fft_relative_windowed = abs(fft_sample_windowed)./abs(fft_reference_windowed);
        phase_relative = unwrap(angle(fft_sample))./unwrap(angle(fft_reference));
        phase_relative_windowed = unwrap(angle(fft_sample_windowed))./unwrap(angle(fft_reference_windowed));
        if zwischenplot
            figure(temp*10+5)
            semilogy(freq_sample,spec_reference_on_air,'k')
            hold on
            semilogy(freq_sample,spec_reference_on_air_windowed,'b')
            
            
            if averaging_in_time
                plot(freq_sample,spec_sample_on_air,'r')
            else
                index = find(freq_sample>ffmax+1,1);
                h = errorbar(freq_sample(1:10:index),spec_sample_on_air(1:10:index),referenced_spec_error(1:10:index),'r')
                h.CapSize = 12;
            end
            
            
            
            
            semilogy(freq_sample,spec_sample_on_air_windowed,'m')
            legend('teflon normalized on air','teflon normalized on air window','Sample normalized on air','Sample normalized on air window')
            xlim([freqmin freqmax])
        end
        %% saving
        endsave=find(freq_sample>4,1,'first');
        fft_relative_windowed_save = fft_relative_windowed(1:endsave);
        fft_relative_save = fft_relative(1:endsave);
        phase_relative_windowed_save= phase_relative_windowed(1:endsave);
        phase_relative_save = phase_relative(1:endsave);
        spectral_intensity_absorbtion_windowed_save = spectral_intensity_absorbtion_windowed(1:endsave);
        spectral_intensity_absorbtion_save = spectral_intensity_absorbtion(1:endsave);
        freq_sample_save = freq_sample(1:endsave);
        if averaging_in_time
        else
            referenced_spec_errorteflon_save = referenced_spec_errorteflon(1:endsave);
            spec_error_save = spec_error(1:endsave);
            save([ordnerspeichern,num2str(temp),'K_spec_error'],'spec_error_save','-ascii');
            save([ordnerspeichern,num2str(temp),'K_referenced_spec_errorteflon'],'referenced_spec_errorteflon_save','-ascii');
        end
        
        
        save([ordnerspeichern,num2str(temp),'K_windowed'],'fft_relative_windowed_save','-ascii');
        save([ordnerspeichern,num2str(temp),'K'],'fft_relative_save','-ascii');
        save([ordnerspeichern,num2str(temp),'K_windowed_phase'],'phase_relative_windowed_save','-ascii');
        save([ordnerspeichern,num2str(temp),'K_phase'],'phase_relative_save','-ascii');
        save([ordnerspeichern,num2str(temp),'K_windowed_spec'],'spectral_intensity_absorbtion_windowed_save','-ascii');
        save([ordnerspeichern,num2str(temp),'K_spec'],'spectral_intensity_absorbtion_save','-ascii');
        save([ordnerspeichern,num2str(temp),'K_freq'],'freq_sample_save','-ascii');
        save([ordnerspeichern,num2str(temp),'zusat'],'probenname','-ascii');
        
        
        probenname
        %    clear A_reference_windowed
        %   clear A_reference
        %   clear freq_sample
        %  clear fenster
    end
    %% Plot commands
    if zwischenplot
        freq_sample = freq_air;
        figure(temp*100+3)
        semilogy(freq_reference,spec_reference,'k')
        hold on
        semilogy(freq_reference,spec_reference_windowed,'b')
        xlim([freqmin freqmax])
        semilogy(freq_sample,spec_sample,'r')
        semilogy(freq_sample,spec_sample_windowed,'m')
        
        figure(temp*100+9)
        plot(freq_air,abs(fft_air))
        hold on
        xlim([freqmin freqmax])
        plot(freq_air,abs(fft_air_windowed),'r')
    end
    %% The following part is only relevant when you would like to use this program to plot previously calculated data
else %% displaying of multiple measurements. Make sure that they have bin evaluated!
    [~,~,~,~,~,~,~,~,ordnerspeichern,~,~,~]= Assigment_sheet_TDS(nummerderreihe,temp);
    if averaging_in_time
        zustatz ='TD';
        ploterrorbar =0;
    else
        zustatz ='FD';
    end
    if real_imag_seperation
        zz ='\';
    else
        zz='\numerical\';
    end
    ordnerspeichern = 'C:\ExperimentalData_backuped\Jens\DLNorleucine_Doped\JN\';
    for lauf =1:length(temprange)
        temp = temprange(lauf);
        
        complexn_windowed = load([ordnerspeichern,num2str(temp),'K_windowed']);
        complexn = load([ordnerspeichern,num2str(temp),'K']);
        freq = load([ordnerspeichern,num2str(temp),'K_freq']);
        specabsorbtion = load([ordnerspeichern,num2str(temp),'K_spec']);
        specabsorbtion_windowed = load([ordnerspeichern,num2str(temp),'K_windowed_spec']);
        phase = load([ordnerspeichern,num2str(temp),'K_phase']);
        phase_windowed = load([ordnerspeichern,num2str(temp),'K_windowed_phase']);
        normalizedon = load([ordnerspeichern,num2str(temp),'normalizedon']);
        
        if averaging_in_time
            figure(42)
        else
            figure(52)
            spec_error =load([ordnerspeichern,num2str(temp),'K_spec_error'],'-ascii');
            referenced_spec_errorteflon = load([ordnerspeichern,num2str(temp),'K_referenced_spec_errorteflon'],'-ascii');
            
            
        end
        ylim([0.3 1.1])
        hold on
        if ploterrorbar
            errorbar(freq(1:50:end),specabsorbtion_windowed(1:50:end),referenced_spec_errorteflon(1:50:end),'color',farbe(lauf,:))
        else
            plot(freq,specabsorbtion_windowed,'color',farbe(lauf,:))
        end
        if averaging_in_time
            figure(421)
        else
            figure(521)
        end
        hold on
        if ploterrorbar
            errorbar(freq,specabsorbtion,referenced_spec_errorteflon,'color',farbe(lauf,:))
        else
            plot(freq,specabsorbtion,'color',farbe(lauf,:))
        end
        if ploterrorbar
            errorbar(freq,specabsorbtion_windowed,referenced_spec_errorteflon,'color',farbe(lauf,:))
        else
        end
        hold off
        if averaging_in_time
            figure(4242)
        else
            figure(5252)
        end
        hold on
        semilogy(freq,specabsorbtion_windowed,'color',farbe(lauf,:));
        if ploterrorbar
            errorbar(freq,specabsorbtion_windowed,referenced_spec_errorteflon,'color',farbe(lauf,:))
        else
        end
        
        
        
        hold off
        if averaging_in_time
            figure(21)
        else
            figure(26)
        end
        hold on
        plot(freq,unwrap(phase_windowed,1),'color',farbe(lauf,:))
        if averaging_in_time
            figure(211)
        else
            figure(266)
        end
        hold on
        plot(freq,unwrap(phase,1),'color',farbe(lauf,:))
        %complexncombined(lauf,:) = complexn;
        
        tempforplot = temprange(lauf);
        
        probenname = load([ordnerspeichern,num2str(temp),'zusat']);
        zusat= probenname;
        lege{lauf}= [zusat,num2str(tempforplot),' K '];
        if fittenn
            ffmin = freqmin;
            ffmax = freqmax;
            nDLN= load([ordnerspeichern,num2str(temp),'refractiveindexsample']);
            nteflon= load([ordnerspeichern,num2str(temp),'refractiveindexteflon']);
            kDLN  = load([ordnerspeichern,num2str(temp),'absorptionindexsample']);
            kteflon = load([ordnerspeichern,num2str(temp),'absorptionindexteflon']);
            nDLN_effectivmedium = load([ordnerspeichern,num2str(temp),'MG_onair']);
            kDLN_effectivmedium  = load([ordnerspeichern,num2str(temp),'MG_onair_absorb']);
            nDLN_effectivemediumT = load([ordnerspeichern,num2str(temp),'MG_onT']);
            kDLN_effectivemediumT  = load([ordnerspeichern,num2str(temp),'MG_onT_absorb']);
            alpha_effektivmediumT = load([ordnerspeichern,num2str(temp),'alphaMG_onT']);
            alpha_effektivmedium = load([ordnerspeichern,num2str(temp),'alphaMG_onT_absorb']);
            frequenzrange = load([ordnerspeichern,num2str(temp),'freqfit']);
            if averaging_in_time
                figure(4444433)
            else
                figure(5555533)
            end
            hold on
            plot(frequenzrange,kDLN_effectivmedium,'x','color',farbe(lauf,:))
            plot(frequenzrange,kDLN_effectivemediumT,'d','color',farbe(lauf,:))
            if averaging_in_time
                figure(444444)
            else
                figure(555555)
            end
            hold on
            plot(frequenzrange,nDLN_effectivmedium,'x','color',farbe(lauf,:))
            plot(frequenzrange,nDLN_effectivemediumT,'d','color',farbe(lauf,:))
            if averaging_in_time
                figure(4343)
                hold on
                
            else
                figure(5353)
                hold on
                
            end
            plot(frequenzrange,alpha_effektivmedium,'color',farbe(lauf,:))
            % plot(frequenzrange,alpha_effektivmediumT,'d','color',farbe(lauf,:))
            legedichhin{(lauf*2)-1}=['ref on air ',zusat,num2str(tempforplot),' K '];
            legedichhin{(lauf*2)}=['ref on teflon ',zusat,num2str(tempforplot),' K '];
            alphaex = alpha_effektivmedium';
            save([ordnerspeichern,'JNexp\',probenname,num2str(tempforplot),' K.m'],'alphaex')
            save([ordnerspeichern,'JNexp\',probenname,num2str(tempforplot),' K.txt'],'alphaex','-ascii')
            freqex = frequenzrange';
            save([ordnerspeichern,'JNexp\',probenname,'freq.m'],'freqex')
            save([ordnerspeichern,'JNexp\',probenname,'freq.txt'],'freqex','-ascii')
            xlim([ffmin ffmax])
        end
    end
    dDLN = load([ordnerspeichern,num2str(temp),'dDLN']);
    if averaging_in_time
        figure(444444)
    else
        figure(555555)
    end
    box on
    legend(legedichhin);
    xlim([ffmin ffmax])
    title(['retrieved realpart of n, ',probenname,'e, for thickness ',num2str(dDLN),' effective medium'])
    if averaging_in_time
        figure(4444433)
    else
        figure(5555533)
    end
    box on
    legend(legedichhin);
    xlim([ffmin ffmax])
    title(['retrieved imagpart of n,',probenname,', for thickness ',num2str(dDLN),' effective medium'])
    if averaging_in_time
        figure(4343)
    else
        figure(5353)
    end
    box on
    legend(num2str(temprange'),'Location','WestOutside');
    xlim([ffmin ffmax])
    ylabel('absorptioncoef. \alpha [cm^{-1}]')
    xlabel('frequency [THz]')
    title(['retrieved absorptioncoefficient, ',probenname,', for thickness ',num2str(dDLN),' effective medium'])
    ymax =1.1;
    ymin =0.01;
    
    if averaging_in_time
        figure(42)
    else
        figure(52)
    end
    
    legend(lege,'Location','SouthWest')
    xlim([freqmin freqmax])
    ylim([ymin ymax])
    xlabel('frequency [THz]')
    ylabel(['relative transmission, normalized on',normalizedon])
    box on
    title([probenname,'normalized on ',normalizedon,', smoothed with gaussian window'])
    if averaging_in_time
        figure(4242)
    else
        figure(5252)
    end
    legend(lege,'Location','SouthWest')
    xlim([freqmin freqmax])
    ylim([0.01 ymax])
    xlabel('frequency [THz]')
    ylabel(['relative transmission, normalized on ',normalizedon])
    box on
    title([probenname,'normalized on',normalizedon,', smoothed with gaussian window'])
    ax=gca;
    ax.YScale = 'log';
    if averaging_in_time
        figure(421)
    else
        figure(521)
    end
    legend(lege,'Location','SouthWest')
    xlim([freqmin freqmax])
    ylim([ymin ymax])
    xlabel('frequency [THz]')
    ylabel(['relative transmission, normalized on',normalizedon])
    title([probenname,' normalized on ',normalizedon])
    box on
    if averaging_in_time
        figure(21)
    else
        figure(26)
    end
    legend(lege,'Location','SouthWest')
    xlim([freqmin freqmax])
    xlabel('frequency [THz]')
    ylabel('phase')
    title([probenname,'normalized on ',normalizedon,', smoothed with gaussian window'])
    box on
    if averaging_in_time
        figure(211)
    else
        figure(266)
    end
    legend(lege)
    xlim([freqmin freqmax])
    xlabel('frequency [THz]')
    ylabel('phase')
    title([probenname,'normalized on ',normalizedon])
    box on
    
end
freqmittelwert = (freq_for_d_start+freq_for_d_end)/2;
if max(temprange)> 10000
    temprange = temprange -10000;
else
    if max(temprange)> 1000
        temprange = temprange -1000;
    end
end

%% And more Plotting ... 
if fittend
    figure(1)
    hold on
    plot(temprange,real(nDLN_all),'k','Marker',symbole,'LineStyle','none')
    plot(temprange,real(nteflon_all),'r','Marker',symbole,'LineStyle','none')
    legend('DLN','Teflon')
    xlabel('Temperature in [K]')
    ylabel(['real refractive index for,',num2str(freqmittelwert),' fitted with Cauchy'])
    title(['used freqrange',num2str(freq_for_d_start),'THz to',num2str(freq_for_d_end),'THz' ])
    figure(3)
    hold on
    plot(temprange,imag(nDLN_all),'k','Marker',symbole,'LineStyle','none')
    plot(temprange,imag(nteflon_all),'r','Marker',symbole,'LineStyle','none')
    legend('DLN','Teflon')
    xlabel('Temperature in [K]')
    ylabel(['imag refractive index for,',num2str(freqmittelwert),' fitted with Cauchy'])
    title(['used freqrange',num2str(freq_for_d_start),'THz to',num2str(freq_for_d_end),'THz' ])
    figure(2)
    hold on
    plot(temprange,dDLN_all,'k','Marker',symbole,'LineStyle','none')
    plot(temprange,dteflon_all,'r','Marker',symbole,'LineStyle','none')
    legend('DLN','Teflon')
    xlabel('Temperature in [K]')
    ylabel('Thickness of the sample, fitted with Cauchy')
    title(['used freqrange',num2str(freq_for_d_start),'THz to',num2str(freq_for_d_end),'THz' ])
    save([ordnerspeichern,'dicke_basedon',num2str(freq_for_d_start),'to',num2str(freq_for_d_end),'.mat'],'dDLN_all','dteflon_all')
end
if evaluatingone
    if mitdergroessespielen
        if fittenn
            nteflon_compl = nteflon+1i*kteflon;
            epsilonh = nteflon_compl.^2;
            nDLN_effectivmedium_compl = nDLN_effectivmedium +1i*kDLN_effectivmedium;
            epsiloni = nDLN_effectivmedium_compl.^2;
            durchmessergesamt = [0 100 200 300 400 700];
            c0 = 2.99792458e8;
            wellenl = (c0./frequenzrange)*1e-6;
            k0 = (2*pi)./(wellenl*1e-6);
            %Ruppin, R. Evaluation of extended Maxwell-Garnett theories Optics Communications, 2000, 182, 273-279
            for rennschweinchen =1:length(durchmessergesamt)
                durchmesser = durchmessergesamt(rennschweinchen);
                sizeparameter = durchmesser./wellenl;
                deltaeffekt = sizeparameter.^2+ ((2/3)*1i).*sizeparameter.^3;
                zaehler = epsiloni.*(1+2*volumenfaktor)+ 2*epsilonh.*(1-volumenfaktor) + (epsilonh-epsiloni).*(1-volumenfaktor).*deltaeffekt;
                nenner =epsiloni.*(1-volumenfaktor)+ epsilonh.*(2+volumenfaktor) + (epsilonh-epsiloni).*(1-volumenfaktor).*deltaeffekt;
                epsilonmean_DMG = epsilonh.*zaehler./nenner;
                n_DMG = sqrt(epsilonmean_DMG);
                alpha_DMG = 2*imag(n_DMG).*k0*1e-2;
                figure(72)
                hold on
                plot(frequenzrange,alpha_DMG,'+','color',farbe(rennschweinchen,:))
                z= (epsiloni./epsilonh -1).*(2/3 *(1-i*sizeparameter).*exp(i*sizeparameter) -1);
                zaehler = 3*epsilonh.*(1-z) + 2*volumenfaktor*(epsiloni-epsilonh);
                nenner = 3*epsilonh.*(1-z) - volumenfaktor*(epsiloni-epsilonh);
                epsilonmean_IMG = epsilonh.*zaehler./nenner;
                n_IMG = sqrt(epsilonmean_IMG);
                alpha_IMG = 2*imag(n_IMG).*k0*1e-2;
                plot(frequenzrange,alpha_IMG,'d','color',farbe(rennschweinchen,:))
                legendemalwider{rennschweinchen*2-1} = ['DMG ',num2str(durchmesser)];
                legendemalwider{rennschweinchen*2} = ['IMG ',num2str(durchmesser)];
            end
            title('calculations to mimic the red-shift of the resonance')
            ylabel('absorbtionskoeffizient [cm^{-1}]')
            xlabel('frequency [THz]')
            legend(legendemalwider)
        end
    end
end


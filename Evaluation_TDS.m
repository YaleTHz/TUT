%% Written by Jens Neu, Yale University
% Labbook: JN p. 91,98,108,123
% Please read the attached readme file.
% All rights reserved (C) Jens Neu 2019, the software can be used according to GPL-3 licence.
% The software is provided as is without any warranty or reliability on my
% site. Please report bugs to jens.neu@yale.edu. Please do not distribute.

% The newest version can be downloaded from thz.yale.edu; There will also
% be a change.log and a FAQ section, as well as known issues.

% Please keep in mind that in the spirit of proper scientific practice,
% used software ought to be cited in any publication in which it is used.
% This software is described in:
% "Exploring the solid state phase transition in DL-norvaline with
% terahertz spectroscopy", by ” by Jens Neu and Coleen T. Nemes and Kevin
% P. Regan and Michael R.C. Williams and Charles A. Schmuttenmaer published
% in Physical Chemistry Chemical Physics, 20, p276-283, 2017 and in
% “Tutorial: An Introduction to Terahertz Time Domain Spectroscopy
% (THz-TDS)” by Jens Neu and Charles A. Schmuttenmaer, published in Journal
% of Applied Physics, 124,22 2018 both publications provide additional
% information on the functioanlity and the underlying physics.

%% Cleaning workspace and setting general commands
clear all % empties workspace, matlabs says I do not need it, but it is needed
offset = true; %substracted an offset from the data, based on the first datapoints
close all %closing all open figures, avoids plotting new/old data in same figure
tic; %starts to clock
%% pointer to data, do edit
farbe = [ 1, 0, 0; 0, 0, 1; 0, 1, 0; 0, 0, 0; 1, 1, 0; 0, 1, 1; 0, 1, 0.5; 1, 0 , 1; 0.25, 0.5,0; 0.5, 0.5 ,0; ...
    1, 0, 0; 0, 0, 1; 0, 1, 0; 0, 0, 0; 1, 1, 0; 0, 1, 1; 0, 1, 0.5; 1, 0 , 1; 0.25, 0.5,0; 0.5, 0.5 ,0; ...
    1, 0, 0; 0, 0, 1; 0, 1, 0; 0, 0, 0; 1, 1, 0; 0, 1, 1; 0, 1, 0.5; 1, 0 , 1; 0.25, 0.5,0; 0.5, 0.5 ,0; ...
    1, 0, 0; 0, 0, 1; 0, 1, 0; 0, 0, 0; 1, 1, 0; 0, 1, 1; 0, 1, 0.5; 1, 0 , 1; 0.25, 0.5,0; 0.5, 0.5 ,0]; %defines a color map
%% less important
exportindex=0; %token that counts the exported lines and ensure that they are not overwritten
exportindex_F3 = exportindex; %if two different kinds of samples, these two indices F3F4 seperate the results
exportindex_F4 = exportindex_F3;
%% program task selection, all boolean,
% keep in mind that 0=false and 1= true
global lit_val_substrate
global lit_val_substrate_value
lit_val_substrate =false; %assumes a constant value for the refractive index of the substrate, if set true
lit_val_substrate_value = 1.99; %3.418; If the previous token is false, this value is ignored!
global reference_on_air
reference_on_air = 1; % if false it uses the reference measurement instead of the air measurement.
evaluatingone = 1; % if false the program only opens prev. results
zwischenplot =1; %makes a shit load of figures
zwischenplott_dfit =1; %makes figures that illustrate the agreement between input and optimized/calculated
zwischenplott_fit =1; %compares clauclated transmission with measured transmission
sigmaplot =0; %plots conductivity
datenfensternfuerfit =1; %etalon truncation on
gaussfenster =1; %etalon truncation with gauss, if false a square function is used
excelexport =0; % Takes a while, I hate excel ...
sigma_einzel =0; % plots each conductivity in an individual figure
smooth =1;  % moving average window, after calculations, this is mainly to make the final output look nicer ....
fittenn =1; % optimization algorithm for refractive index
druiden =0; % requires fitten; uses drude or drude-smith models
lattice =0; % if a non-conducting sample is also measured, this sample can provide an "static" offset, rarely used ... some examples would be doped silicon and SI-silicon as additional sample...
fittend =1; % fitts the refractive index of the substrate with a model (debye, cauchy, or Sellmeier), needs fitten! It will also fit the thickness of the sample!
fittend_sample =0 ; %fitts the refractive index of the sample with a model (debye, cauchy, or Sellmeier), needs fitten! It will also fit the thickness of the sample!
fittend_no_model =0; % calculates the optimal thickness by minimizing the wavyness of n,k (takes for ever..., I am nut sure if properly implemented!!!)
load_substrat =0; %load the substrate paramater from a prev. generated file. It will pick the closest freq point, not average, so ensure proper data density in the load data
tapecell =0; % supports a sample sandwiched between two substrates, as it would be in a "tapecell"
fit_d_on_reduced_data =1; %limits the frequency range used for fitting/calculations, this is used if the high/low frequency range of your measurement are more noisy and you do not want them in the fit, but still plotted alongside the fit
tempinvert =1; % boolean, flips the time axis if needed
if tapecell % using the reference instead of air is not yet implemented for tape cell, therefore this if/case overwrites the previous selection....
    reference_on_air = 1;
    index2 =3; %index2 can be 1 or 2 if you want to split between two series of samples. If you do not need it, set it on a value larger than 2
else
    index2 =4;
end
global drudesmith
global LDM
drudesmith =0;% requires druiden, uses drude-smith instead of drude
LDM =0;% requires druiden, overwrites drudesmith! Uses a localiced plasmon model, not implemented yet!!!
%% choose measurement, These measurements must be defined in the addigment file
% this program is a two piece program. You need to define each "temp" in a
% seperated file (see readme.txt)
temprange = [3301];
%vector containing all pre-defined measurement IDs
%temprange = 292;

%% control variable for fit, do not edit, most of the alternative functions are not correct and therefore haven't been used
global freq_for_d_start
global freq_for_d_end
freq_for_d_start = 0.35; % frequency range used for thickness fit
freq_for_d_end = 2;% frequency range used for thickness fit
global teflon
global frequenzschritt % frequency stepsize for pointwise calculation of refractive index
frequenzschritt =50; % in ghz
global frequenzmitteln %averaging in frequency domain.
frequenzmitteln =1; %false uses the frequency point at each frequency step, true averages the data inbetween
teflon =1; %this variable is used to pick the correct equations, whether substrate or sample is evaluated. Here it is only initialized
%% defines the frequency range, for refractive index calculation, plotting, and model-fitting
freqmin = 0.35; % frequency range used for plot
freqmax = 3.0;
freqmin_thickness = 0.4; % frequency range used for thickness fit
freqmax_thickness = 2.5;
fanfang =0.25; % and another freq. range for calculations
fende = freqmax;
epsilon_0 = 8.85419e-12; %physical constant, only change this if the speed of light has change and the universe is messed up. If you believe in Maxwell eq. keep this constant
global dpila
global dteflon
dteflon =1180;%1180; %SiO2/substrate thickness in mue, historical reasons.1.18E3, can be overwritten in the other file!
global sellerie
global debyemodel
global cauchy
global lorentzian
sellerie =0; debyemodel =0; cauchy =0; lorentzian =1;% one of these must be true for fitting the thickness and refractive index simultaneously
%are overwritten in the call of the fit program
global mmmmm % one of many control variables allowing data exchange between function and main prog. This ones stores the phase-branch in case 2pi jumps occur.
%% Here the program starts
weite =6.5; % width of time window to supress etalons
if evaluatingone
    for lauflauf =1:length(temprange) %loops over all elements in the temprange, this allows you to process large series
        farbeglobal = hsv(length(temprange)); %creates another color map
        temp = temprange(lauflauf); % defines the current measurement
        %% some stuff for the legend and plot
        corr=0;
        lege{lauflauf*2} = num2str(temp-corr);
        lege{lauflauf*2-1} = 'ref';
        lale{lauflauf} = num2str(temp-corr);
        LW =2; %defines the linewidth of some plots, must be a number in pts
        %% assignmend of measured data, CHANGE THIS!!!
        runID = 'F5,F7 pos0'; % (string) can be used to seperate different series, in the same assigment file
        nummerderreihe =1; % (integer) can be used to seperate different series, in the same assigment file
        nummerderprobe =1; % (integer) can be used to seperate different series, in the same assigment file
        [temprangeall,monat_messung,tag_messung,name_messung_air,name_messung_reference,name_messung_sample,...
            ordner,measurementtype,ordnerspeichern,nummer_messung_sample,nummer_messung_reference,nummer_messung_air,probenname]...
            = Assigment_sheet_TDS(nummerderreihe,temp); %retrieves the parameter from the assigment file
        %% transfer to fit program
        dorig_teflon = dteflon;
        global freq_air
        global freq_sample
        global freq_reference
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
        %% opening the data
        %% assigning data
        [t_air, A_air] = Import_for_timfiles(name_messung_air);
        [t_reference, A_reference] = Import_for_timfiles(name_messung_reference);
        [t_sample, A_sample] = Import_for_timfiles(name_messung_sample);
        tchange = 70.5; % this is used for a binary window on the data. it is overwritten by the Gausswindow command in the main settings part
        %% plot input time data
        figure(989898)
        hold on
        plot(t_air,A_air,'c')
        plot(t_reference,A_reference,'b')
        plot(t_sample,A_sample,'m')
        %% stores the different names in one single vector/array to allow usage for labels and legends
        name_air{lauflauf} = name_messung_air;
        name_ref{lauflauf} = name_messung_reference;
        name_sample{lauflauf} = name_messung_sample;
        %% if gaussian window is NOT used, this defines Stepfunction-windows
        for lauf=1:length(t_sample)
            if t_sample(lauf) > tchange
                sprungfenster_sample(lauf,1) =0;
            else
                sprungfenster_sample(lauf,1) =1;
            end
        end
        for lauf=1:length(t_air)
            if t_air(lauf) > tchange
                sprungfenster_air(lauf,1) =0;
            else
                sprungfenster_air(lauf,1) =1;
            end
        end
        [maxi index] = max(A_sample);
        deltat =t_sample(index);
        %% windows data
        if gaussfenster
            fenster_sample = exp(-((t_sample-deltat)/weite).^2);
        else
            fenster_sample = sprungfenster_sample;
        end
        A_sample_windowed = A_sample.*fenster_sample;
        [maxi index] = max(A_reference);
        deltat =t_reference(index);
        if gaussfenster
            fenster = exp(-((t_reference-deltat)/weite).^2);
        else
            fenster = sprungfenster_sample;
        end
        [maxi index] = max(A_air);
        deltat =t_sample(index);
        if gaussfenster
            fenster_air = exp(-((t_air-deltat)/weite).^2);
        else
            fenster_air = sprungfenster_air;
        end
        
        A_reference_windowed = A_reference.*fenster;
        %% some FFT
        NNT = 2^18;% better 15 or 16, Zero padding!
        freq_reference = faxis(t_reference,NNT); %faxis is a small function that creates THz axis if the input data is in ps unit
        spec_reference = abs(fft(A_reference,NNT));
        fft_reference = fft(A_reference,NNT);
        spec_reference_windowed = abs(fft(A_reference_windowed,NNT));
        fft_reference_windowed = fft(A_reference_windowed,NNT);
        freq_sample = faxis(t_sample,NNT);
        spec_sample = abs(fft(A_sample,NNT));
        spec_sample_windowed = abs(fft(A_sample_windowed,NNT));
        fft_sample = fft(A_sample,NNT);
        fft_sample_windowed = fft(A_sample_windowed,NNT);
        freq_air = faxis(t_air,NNT);
        fft_air = fft(A_air,NNT);
        A_air_windowed = A_air.*fenster_air;
        fft_air_windowed = fft(A_air_windowed,NNT);
        spec_reference_on_air = abs(fft_reference)./abs(fft_air);
        spec_reference_on_air_windowed = abs(fft_reference_windowed)./abs(fft_air_windowed);
        spec_sample_on_air = abs(fft_sample)./abs(fft_air);
        spec_sample_on_air_windowed = abs(fft_sample_windowed)./abs(fft_air_windowed);
        S21_pila = fft_sample_windowed./fft_air_windowed;
        S21_SiO2 = fft_reference_windowed./fft_air_windowed;
        freq = freq_air;
        dpila_max =100;
        dpila_min =0;
        %% Approximation of n from timedomain, used as starting point for optimization
        ddpila = dpila*1e-6;  % in meter from here on
        ddsio2 =dteflon*1e-6; % in meter from here on
        c=physconst('LightSpeed');
        [maxi index] = max(A_air);
        tmaxair =t_air(index);
        [maxi index] = max(A_reference);
        tmaxref =t_reference(index);
        [maxi index] = max(A_sample);
        tmaxpila =t_sample(index);
        deltasio2 = (tmaxref-tmaxair)*1e-12;
        nmean_sio2 = ((c*deltasio2)/ddsio2)+1;
        deltapila = (tmaxpila-tmaxref)*1e-12;
        nmean_pila = ((c*deltapila)/ddpila)+1;
        %% preparing data for optimization: Discretize data, ensure correct phase, assign to global var.
        global frequenzrange
        frequenzrange = fanfang:(frequenzschritt/1000):fende;
        omega = 2*pi*frequenzrange*1e12;
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
        if tempinvert
            Esample = conj(Esample);
            Eair = conj(Eair);
            Ereference = conj(Ereference);
        end
        global S21teflon
        S21teflon = Ereference./Eair;
        global S21DLN
        S21DLN = Esample./Eair;
        global S21prime
        S21prime = Esample./Ereference;
        global S21teflonphase
        global S21DLNphase
        global S21DLN_on_Teflon_phase
        global S21DLN_on_Teflon
        global S21target
        S21DLN_on_Teflon = Esample./Ereference;
        S21DLN_on_Teflon_phase = unwrap(angle(S21DLN_on_Teflon));
        S21teflonphase = unwrap(angle(S21teflon));
        S21DLNphase = unwrap(angle(S21DLN));
        global nteflon
        global kteflon
        %% calling the fit functions
        if fittenn
            global specnorm
            specnorm = max(abs(Eair)); %originally used to normalize the optimization goals, I am fairly certain not in use anymore, however as these are global variables I am not brave enough to delete them ..
            global phasenorm
            phasenorm = max(unwrap(angle((Eair))));
            global frequenzpunkt
            global noisecancelling
            noisecancelling = 1; % averaging in GHz in freq. domain, if set 0 no averaging. Not to be confused with smoothing, which is applied after the fit
            % seed for n and k
            nteflon(1) = 2;%3.418;
            kteflon(1) = 0;
            global S21_fuer_dopti
            global S21_fuer_dopti_teflon
            for fitschleife =2:length(frequenzrange)
                frequenzpunkt = frequenzrange(fitschleife);
                % refractive index of mixture and reference
                teflon =1;
                if frequenzmitteln %averages signal data in the frequency domain
                    range = frequenzschritt/2000;% /1000 for GHz->THz, /2 for half +/-
                    f = frequenzpunkt;
                    argh = f-range;
                    [indexair_lb mini] = find(freq_air>argh,1,'first');
                    [indexsample_lb mini] = find(freq_sample>argh,1,'first');
                    [indexreference_lb mini] = find(freq_reference>argh,1,'first');
                    argh = f+range;
                    [indexair_ub mini] = find(freq_air>argh,1,'first');
                    [indexsample_ub mini] = find(freq_sample>argh,1,'first');
                    [indexreference_ub mini] = find(freq_reference>argh,1,'first');
                    Slocal  = mean(S21DLN(indexreference_lb:indexreference_ub));
                    Slocalteflon = mean(S21teflon(indexreference_lb:indexreference_ub));
                else
                    [indexair mini] = find(freq_air>f,1,'first');
                    [indexsample mini] = find(freq_sample>f,1,'first');
                    [indexreference mini] = find(freq_reference>frequenzpunkt,1,'first');
                    Slocal = S21DLN(indexsample);
                    Slocalteflon = S21teflon(indexsample);
                end
                S21_fuer_dopti(fitschleife) = Slocal;
                S21_fuer_dopti_teflon(fitschleife) = Slocalteflon;
                
            end
            
            
        end
        global freq_dfit
        if fit_d_on_reduced_data
            freq_dfit=frequenzrange;
        else
            startindex= find(freq_reference>freqmin_thickness,1,'first');
            endindex = find(freq_reference>freqmax_thickness,1,'first');
            freq_dfit = freq_reference(startindex:endindex);
            S21_fuer_dopti_teflon = S21teflon(startindex:endindex);
        end
        %% fitting thickness
        if fittend
            clear x %because x has to clear from time to time ....
            teflon =1; %ensures that the equations for the substrate are used, if your substrate is Silicon or glass, keep it still on true, the name teflon is historical
            x(1) = dteflon; %starting value
            %% These values depend on the chosen model, double check which to use.
            %% asigns input values from easier to read physical quanteties to numerically more stable values around -1 to 1...
            if debyemodel
                x(2) = 3; %debye model, epsilon inf
                x(3) = 1; %debye model, epsilon res1
                x(4) = 7; %debye model, tau1
                x(5) = 0.4; %debye model epsilon res2
                x(6) = 0.12; %debye model tau2
                x(7) = 0.2; %debye model epsilon res2
                x(8) = 0.12; %debye model tau2
                x(9) = 0; %debye model epsilon res2
                x(10) = 0; %debye model tau2
                lb = x*0.0; ub = x*1000000;
            end
            if cauchy
                x(2) = 1.9; % A, constant value part
                x(3) = 0.0001; %B first resonant term
                x(4) = 0.00001; %C lambda4 term
                x(5) = 0.000; % alpha, absorption strnght
                x(6) = 1; % beta exponential strength
                x(7) = 100; % gamma absorbing resonance, in µm
                lb = x*0.0001; ub = x*10000;
            end
            if sellerie
                x(2)= 1.9;%first resonant term
                x(3)= 300; %resonance wavelength in µm
                x(4)= 2; %second resonance term
                x(5)= 600; %second wavelength
                x(6) = 0;% absorbtion, constant.
                lb = x*0.0001; ub = x*10000;
            end
            if lorentzian
                x(2) = 1.84; % epsilon inf
                x(3) = 6.1; % omega_resonance  *2pi in THz
                x(4) = 0.6; % Omega_resonance  *2pi in THz
                x(5) = 0.2; % gamma, damping
            end
            lb(1) = dorig_teflon*0.75; ub(1) = dorig_teflon*1.25;
            
            options = optimoptions('fmincon');
            options = optimoptions(options,'Display', 'off');
            options = optimoptions(options,'FunValCheck', 'off');
            [x,fval,exitflag,output,lambda,grad,hessian] =...
                fmincon(@retrieval_of_d_mark3,x,[],[],[],[],lb,ub,[],options); %calls one of the optimization algorithms
            dteflon = x(1);
            %% reasigns optimized values to easier to read physical quanteties, and replots result
            if debyemodel
                eps_inf =x(2);
                eps_res1 =x(3);
                tau1=x(4)*1E-12;
                eps_res2 =x(5);
                tau2=x(6)*1E-12;
                eps_res3 =x(7);
                tau3=x(8)*1E-12;
                eps_res4 =x(9);
                tau4=x(10)*1E-12;
                complexomega = 2*pi*1i* frequenzrange*1E12;
                epsilon_deb = eps_inf +(eps_res1./(1+complexomega.*tau1))  +(eps_res2./(1+complexomega*tau2))+(eps_res3./(1+complexomega*tau3))+(eps_res4./(1+complexomega*tau4));
                n_model_sio = sqrt(epsilon_deb);
            end
            if cauchy
                A =x(2); B =x(3); C =x(4); alpha =x(5); beta= x(6); gamma = x(7);
                c0 = 2.99792458e8;
                wellenl = (c0./frequenzrange)*1e-6; %wavelength in micrometer
                ll = wellenl.*wellenl;
                llll = ll.*ll;
                n_model_sio = A + B./ll+C./llll + 1i.*(alpha.*exp(12400.*(beta.*(1./wellenl- 1./gamma))));
            end
            if sellerie
                B1 = x(2);C1 = x(3); B2 = x(4); C2 = x(5); absorber = x(6);
                c0 = 2.99792458e8;
                wellenl = (c0./frequenzrange)*1e-6; %wavelength in micrometer
                ll = wellenl.*wellenl;
                llll = ll.*ll;
                n_model_sio = sqrt(1+ B1*ll./(ll-C1) + B2*ll./(ll-C1)) + absorber*1i;
            end
            if lorentzian
                eps_inf = x(2);
                omega_1 = x(3)*2*pi; %all in THz
                Omega_1 = x(4) *2*pi;
                gamma_1 = x(5);
                omega = frequenzrange*2*pi;
                zaehler = (omega_1.^2) - (omega.^2) - 1i*gamma_1*omega;
                epsilon = eps_inf + (Omega_1.^2)./(zaehler);
                n_model_sio = sqrt(epsilon);
            end
            x_siO =x;
            %% plots a comparison of input data and optimization output. If both look identical the optimization might be trustworthy...
            if zwischenplott_dfit
                global freq_dfit_result
                global S21_input
                global S21_output
                figure(140+lauflauf)
                hold on
                plot(freq_dfit_result,abs(S21_input),'r')
                plot(freq_dfit_result,abs(S21_output),'k')
                plot(frequenzrange,abs(S21_fuer_dopti_teflon),'b:')
                xlim([freq_for_d_start freq_for_d_end])
                figure(170+lauflauf)
                hold on
                plot(freq_dfit_result,unwrap(angle(S21_input)),'r')
                plot(freq_dfit_result,unwrap(angle(S21_output)),'k')
                plot(frequenzrange,unwrap(angle(S21_fuer_dopti_teflon)),'b:')
                xlim([freq_for_d_start freq_for_d_end])
            end
        end
        %% fitting the sample
        if fittenn
            if load_substrat %the program can use previous substrate data.
                nteflon_all = load('C:\ExperimentalData_backuped\brian\Brian\Conductive MOF\TDS_results\fixed_d\tapecell\mean_nmean','-ascii');
                kteflon_all = load('C:\ExperimentalData_backuped\brian\Brian\Conductive MOF\TDS_results\fixed_d\tapecell\mean_kmean','-ascii');
                nteflon_all_freq = load('C:\ExperimentalData_backuped\brian\Brian\Conductive MOF\TDS_results\fixed_d\tapecell\292K_frequenzrange_forallfits.txt','-ascii');
                for fitschleife =2:length(frequenzrange)
                    frequenzpunkt = frequenzrange(fitschleife);
                    [~,load_index] = min(abs(nteflon_all_freq-frequenzpunkt));
                    nteflon(fitschleife)=nteflon_all(load_index);
                    kteflon(fitschleife)=kteflon_all(load_index);
                end
            else  %or calculate the substrate parameters based on the reference and air measurements.
                global S21target_teflon
                nteflon(1) = 1.4;
                kteflon(1) = 0;
                for fitschleife =2:length(frequenzrange)
                    frequenzpunkt = frequenzrange(fitschleife);
                    teflon =1;
                    if tapecell
                        x(1) = nteflon(fitschleife-1);x(2) = kteflon(fitschleife-1);
                        lb = [1 0]; ub = [2 2];
                        options = optimoptions('fmincon');
                        options = optimoptions(options,'Display', 'off');
                        options = optimoptions(options,'FunValCheck', 'off');
                        [x,fval,exitflag,output,lambda,grad,hessian] =...
                            fmincon(@retrieval_of_n_airref_tapecell,x,[],[],[],[],lb,ub,[],options);
                        optimization(fitschleife)=S21target_teflon;
                        nteflon(fitschleife) = x(1);kteflon(fitschleife) = x(2);
                        mtest_teflon(fitschleife) =mmmmm;
                    else
                        [nteflon(fitschleife),kteflon(fitschleife)] = retrieval_of_n_airref_nooptimization(nteflon(fitschleife-1),kteflon(fitschleife-1));
                        S21teflon_processed(fitschleife) = abs(S21target);
                        
                    end
                end
                if tapecell
                    figure(430+lauflauf)
                    hold on
                    plot(frequenzrange,abs(optimization),'r')
                    plot(freq_reference,abs(S21teflon),'k')
                    xlim([0.3 2.5])
                    figure(460+lauflauf)
                    hold on
                    plot(frequenzrange,angle(optimization),'r')
                    plot(freq_reference,unwrap(angle(S21teflon)),'k')
                    plot(frequenzrange,mtest_teflon,'b')
                    xlim([0.3 2.5])
                end
                nteflon(1) = nteflon(2);
                kteflon(1) = kteflon(2);
            end
        end
        if fittend & fittend_sample %basically the same as fittend, but now the sample and not the substrate is modelled and fitted
            teflon =0; % this switch to sample calculations!
            x(1) = dpila; %starting value
            %% These values depend on the chosen model, double check which to use.
            if debyemodel
                x(2) = 3; %debye model, epsilon inf
                x(3) = 1; %debye model, epsilon res1
                x(4) = 7; %debye model, tau1
                x(5) = 0.4; %debye model epsilon res2
                x(6) = 0.12; %debye model tau2
                x(7) = 0.2; %debye model epsilon res2
                x(8) = 0.12; %debye model tau2
                x(9) = 0; %debye model epsilon res2
                x(10) = 0; %debye model tau2
                lb = x*0.0; ub = x*1000000;
            end
            if cauchy
                x(2) = 1.9; % A, constant value part
                x(3) = 0.0001; %B first resonant term
                x(4) = 0.00001; %C lambda4 term
                x(5) = 0.000; % alpha, absorption strnght
                x(6) = 1; % beta exponential strength
                x(7) = 100; % gamma absorbing resonance, in µm
                lb = x*0.0001; ub = x*10000;
            end
            if sellerie
                x(2)= 1.9;%first resonant term
                x(3)= 300; %resonance wavelength in µm
                x(4)= 2; %second resonance term
                x(5)= 600; %second wavelength
                x(6) = 0;% absorbtion, constant.
                lb = x*0.0001; ub = x*10000;
            end
            if lorentzian
                x(2) = 1.84; % epsilon inf
                x(3) = 6.1; % omega_resonance  *2pi in THz
                x(4) = 0.6; % Omega_resonance  *2pi in THz
                x(5) = 0.2; % gamma, damping
            end
            lb = x*0.025; ub = x*40;
            lb(1) = dpila_min; ub(1) =dpila_max;x(1)=dpila;
            options = optimoptions('fmincon');
            options = optimoptions(options,'Display', 'off');
            options = optimoptions(options,'FunValCheck', 'off');
            [x,fval,exitflag,output,lambda,grad,hessian] =...
                fmincon(@retrieval_of_d_mark3,x,[],[],[],[],lb,ub,[],options);
            dpila = x(1);
            
            if debyemodel
                eps_inf =x(2);
                eps_res1 =x(3);
                tau1=x(4)*1E-12;
                eps_res2 =x(5);
                tau2=x(6)*1E-12;
                eps_res3 =x(7);
                tau3=x(8)*1E-12;
                eps_res4 =x(9);
                tau4=x(10)*1E-12;
                complexomega = 2*pi*1i* frequenzrange*1E12;
                epsilon_deb = eps_inf +(eps_res1./(1+complexomega.*tau1))  +(eps_res2./(1+complexomega*tau2))+(eps_res3./(1+complexomega*tau3))+(eps_res4./(1+complexomega*tau4));
                n_debye = sqrt(epsilon_deb);
            end
            if cauchy
                A =x(2); B =x(3); C =x(4); alpha =x(5); beta= x(6); gamma = x(7);
                c0 = 2.99792458e8;
                wellenl = (c0./frequenzrange)*1e-6; %wavelength in micrometer
                ll = wellenl.*wellenl;
                llll = ll.*ll;
                n_debye = A + B./ll+C./llll + 1i.*(alpha.*exp(12400.*(beta.*(1./wellenl- 1./gamma))));
            end
            if sellerie
                B1 = x(2);C1 = x(3); B2 = x(4); C2 = x(5); absorber = x(6);
                c0 = 2.99792458e8;
                wellenl = (c0./frequenzrange)*1e-6; %wavelength in micrometer
                ll = wellenl.*wellenl;
                llll = ll.*ll;
                n_debye = sqrt(1+ B1*ll./(ll-C1) + B2*ll./(ll-C1)) + absorber*1i;
            end
            if lorentzian
                eps_inf = x(2);
                omega_1 = x(3)*2*pi; %all in THz
                Omega_1 = x(4) *2*pi;
                gamma_1 = x(5);
                omega = frequenzrange*2*pi;
                zaehler = (omega_1.^2) - (omega.^2) - 1i*gamma_1*omega;
                epsilon = eps_inf + (Omega_1.^2)./(zaehler);
                n_debye = sqrt(epsilon);
            end
        end
        
        %% plotting, a lot ...
        
        if zwischenplot
            figure(temp*10+1)
            hold on
            plot(t_air,A_air,'c')
            plot(t_reference,A_reference,'k')
            plot(t_reference,A_reference_windowed,'b')
            plot(t_sample,A_sample,'r')
            plot(t_sample,A_sample_windowed,'m')
            plot(t_sample,fenster_sample*0.5e-2,'g')
            legend('A_air','A_reference','A_reference_windowed','A_sample','A_sample_windowed','fenster')
            title('Time signal')
            figure(temp*10+2)
            hold on
            plot(freq_reference,spec_reference,'k')
            plot(freq_reference,spec_reference_windowed,'b')
            xlim([freqmin freqmax])
            plot(freq_sample,spec_sample,'r')
            plot(freq_sample,spec_sample_windowed,'m')
            plot(freq_air,abs(fft_air),'c')
            plot(freq_air,abs(fft_air_windowed),'g')
            legend('spec_ref','spec_ref_window','spec_sample','spec_sample_window','air','air windowed')
            title('Linear Amplitude plot')
        end
        spectral_intensity_absorbtion = spec_sample./spec_reference;
        spectral_intensity_absorbtion_windowed = spec_sample_windowed./spec_reference_windowed;
        
        if zwischenplot
            figure(temp*10+4)
            plot(freq_sample,spectral_intensity_absorbtion,'k')
            hold on
            plot(freq_sample,spectral_intensity_absorbtion_windowed,'b')
            xlim([freqmin freqmax])
            ylim([0.3 1.1])
            legend('spectral_intensity_absorbtion','spectral_intensity_absorbtion_windowed')
            title('sample/reference')
        end
        fft_relative = abs(fft_sample)./abs(fft_reference);
        fft_relative_windowed = abs(fft_sample_windowed)./abs(fft_reference_windowed);
        phase_relative = unwrap(angle(fft_sample))./unwrap(angle(fft_reference));
        phase_relative_windowed = unwrap(angle(fft_sample_windowed))./unwrap(angle(fft_reference_windowed));
        
        if zwischenplot
            figure(temp*10+5)
            %semilogy(freq_sample,spec_reference_on_air,'k')
            semilogy(freq_sample,spec_reference_on_air_windowed,'b:')
            hold on
            %semilogy(freq_sample,spec_sample_on_air,'r')
            semilogy(freq_sample,spec_sample_on_air_windowed,'m:')
            xlim([freqmin freqmax])
            legend('spec_reference_on_air','spec_reference_on_air_windowed','spec_sample_on_air','spec_sample_on_air_windowed')
            title('Log amplitude plot')
        end
        if zwischenplot
            figure(temp*10+3)
            semilogy(freq_reference,spec_reference,'k')
            hold on
            semilogy(freq_reference,spec_reference_windowed,'b')
            xlim([freqmin freqmax])
            semilogy(freq_sample,spec_sample,'r')
            semilogy(freq_sample,spec_sample_windowed,'m')
            semilogy(freq_air,abs(fft_air_windowed),'c')
            title('Log sample and reference')
            figure(temp*10+9)
            plot(freq_air,abs(fft_air))
            hold on
            xlim([freqmin freqmax])
            plot(freq_air,abs(fft_air_windowed),'r')
            title('Spectral Air')
            
        end
        tempsaveing = mod(temp,1000);
        %% saving, creates the subfolders, however you need to create the main folder by hand ...
        if fittend
            probenname = ['fitted_d\',probenname];
        else
            probenname = ['fixed_d\',probenname];
        end
        if (exist([ordnerspeichern,probenname,'\'])) ==7
        else
            mkdir([ordnerspeichern,probenname,'\'])
        end
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_windowed'],'fft_relative_windowed','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K'],'fft_relative','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_windowed_phase'],'phase_relative_windowed','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_phase'],'phase_relative','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_windowed_spec'],'spectral_intensity_absorbtion_windowed','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_spec'],'spectral_intensity_absorbtion','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_freq'],'freq_sample','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_t_sample.txt'],'t_sample','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_A_sample.txt'],'A_sample','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_SiO_freq.txt'],'freq_reference','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_SiO_spec.txt'],'spec_reference','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_t_air.txt'],'t_air','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_A_air.txt'],'A_air_windowed','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_t_ref.txt'],'t_reference','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_A_ref.txt'],'A_reference_windowed','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'thickness.txt'],'dpila','-ascii');
        %transmission_sample = spectral_intensity_absorbtion_windowed./spec_reference;
        %transmission_SiO = spectral_intensity_absorbtion_windowed./spec_reference;
        speciair= abs(Eair);
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_spec_air.txt'],'speciair','-ascii');
        
        
        % else %% displaying of multiple measurements. Make sure that they have bin evaluated!
        %     close all
        %     for lauf =1:length(temprange)
        %         temp = temprange(lauf);
        %         complexn_windowed = load([ordnerspeichern,num2str(temp),'K_windowed']);
        %         complexn = load([ordnerspeichern,num2str(temp),'K']);
        %         freq = load([ordnerspeichern,num2str(temp),'K_freq']);
        %         specabsorbtion = load([ordnerspeichern,num2str(temp),'K_spec']);
        %         specabsorbtion_windowed = load([ordnerspeichern,num2str(temp),'K_windowed_spec']);
        %         phase = load([ordnerspeichern,num2str(temp),'K_phase']);
        %         phase_windowed = load([ordnerspeichern,num2str(temp),'K_windowed_phase']);
        %         figure(42)
        %         ylim([0.3 1.1])
        %         hold on
        %         plot(freq,specabsorbtion_windowed,'color',farbe(lauf,:))
        %         figure(421)
        %         hold on
        %         plot(freq,specabsorbtion,'color',farbe(lauf,:))
        %         ylim([0.3 1.1])
        %         figure(21)
        %         hold on
        %         plot(freq,unwrap(phase_windowed,1),'color',farbe(lauf,:))
        %         figure(211)
        %         hold on
        %         plot(freq,unwrap(phase,1),'color',farbe(lauf,:))
        %         %complexncombined(lauf,:) = complexn;
        %         lege{lauf}= [num2str(temprange(lauf)),' K '];
        %     end
        
        
        
        %% more plotting ....
        if zwischenplot
            figure(42)
            %   legend(lege,'Location','SouthWest')
            xlim([freqmin freqmax])
            xlabel('frequency [THz]')
            ylabel('relative transmission, normalized on Teflon')
            box on
            title('normalized on teflon, smoothed with gaussian window')
        end
        if index2 ==1;
            symsymsym =':';
        else
            symsymsym = '*';
        end
        S21_phase_pila = unwrap(angle(S21_pila));
        S21_abs_pila = abs(S21_pila);
        S21_phase_SiO2 = unwrap(angle(S21_SiO2));
        S21_abs_SiO2 = abs(S21_SiO2);
        c=physconst('LightSpeed');
        index3 =0;
        nu = frequenzrange;
        omega = 2*pi*nu*1e12;
        k_0 = omega/c;
        symsym ='-';
        sym ='-';
        if mod(temp,25) ==0
            sym =':';
            symsym ='+';
        end
        if mod(temp,25) ==1
            if index2 ==1
                sym =':';
            else
                sym ='-';
            end
            symsym ='x';
        end
        if mod(temp,25) ==2
            sym ='-';
            symsym ='d';
        end
        temp_farbe = mod(temp,1000);
        if temp_farbe>470 & temp_farbe < 600
            farbelocal = [(temp_farbe-470)/130 0 1-((temp_farbe-470)/130)];
        else
            farbelocal =farbeglobal(lauflauf,:);
        end
        %% and even more plotting
        if zwischenplot
            figure(42)
            plot(freq,S21_phase_SiO2,'k')
            hold on
            plot(freq,S21_phase_pila,sym,'color',farbelocal,'LineWidth',LW)
            xlim([0.4 2.5])
            figure(99)
            hold on
            plot(freq,S21_phase_pila-S21_phase_SiO2,sym,'color',farbelocal,'LineWidth',LW)
            xlim([freqmin freqmax])
            legend('Quartz','PilA')
            title('Phase normalized on air')
            figure(43)
            plot(freq,S21_abs_SiO2,'k')
            hold on
            plot(freq,S21_abs_pila,symsymsym,'color',farbelocal,'LineWidth',LW)
            xlim([freqmin freqmax])
            legend('Quartz','PilA')
            title('Spectraltrans normalized on air')
            ylabel('Transmission')
            xlabel('frequency [THz]')
        end
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_S21_abs_SiO2.txt'],'S21_abs_SiO2','-ascii');
        save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_S21_abs_pila.txt'],'S21_abs_pila','-ascii');
        
        figure(44)
        hold on
        plot(freq,S21_abs_pila./S21_abs_SiO2,symsymsym,'color',farbelocal,'LineWidth',LW)
        xlim([freqmin freqmax])
        legend('Quartz','PilA')
        title('Spectrum normalized on SiO2')
        
        figure(45)
        hold on
        plot(freq,S21_abs_SiO2,sym,'color',farbelocal,'LineWidth',LW)
        xlim([freqmin freqmax])
        title('Spectrum of SiO2 versus temp')
        
        tempi = mod(temp,1000);
        legend_for_43{lauflauf*2} = [probenname,'at',num2str(tempi)];
        legend_for_43{lauflauf*2-1} = ['ref'];
        legend_for_44{lauflauf}= [probenname,'at',num2str(tempi)];
        %% and so much more plotting ....
        if fittenn
            figure(333)
            hold on
            plot(frequenzrange,kteflon,'color',farbelocal)
            if fittend
                plot(frequenzrange,imag(n_model_sio),':','color',farbelocal)
            end
            box on
            title('imag(n)')
            xlabel('frequency [THz]')
            box on
            figure(666)
            hold on
            plot(frequenzrange,nteflon,'color',farbelocal)
            if fittend
                plot(frequenzrange,real(n_model_sio),':','color',farbelocal)
                dteflonnn(lauflauf) = dteflon;
            end
            title('real(n)')
            xlabel('frequency [THz]')
            %% some more exporting ...
            save([ordnerspeichern,'\export\',num2str(temp),'tesa_n'],'nteflon','-ascii')
            save([ordnerspeichern,'\export\',num2str(temp),'tesa_k'],'kteflon','-ascii')
            save([ordnerspeichern,'\export\',num2str(temp),'tesa_freq'],'frequenzrange','-ascii')
            %% and plotting again
            box on
            alphateflon = (2.*kteflon.*(frequenzrange.*2.*pi*1e12)./c)./100;
            figure(999)
            hold on
            plot(frequenzrange,alphateflon,'color',farbelocal)
            title('absorptioncoefficient')
            xlabel('frequency [THz]')
            ylabel('absorption /cm')
            box on
            specnorm = max(abs(Eair));
            phasenorm = max(unwrap(angle((Eair))));
            noisecancelling = 1; % averaging in GHz in freq. domain, if set 0 no averaging.
            % seed for n and k
            nDLN(1)  = nmean_pila; %might work better.....
            
            if index2==1
                kDLN(1) =0.5;
            else
                kDLN(1) =5;
            end
            %% calculating the sample refractive index
            for fitschleife =2:length(frequenzrange)
                
                global n_substrat_loc
                n_substrat_loc=nteflon(fitschleife)-1i*kteflon(fitschleife);
                frequenzpunkt = frequenzrange(fitschleife);
                % refractive index of mixture and reference
                teflon =0;
                x(1) = nDLN(fitschleife-1)+0;
                x(2) =kDLN(fitschleife-1)+0;
                %teflon_thickness, refrective index no dispersion, absorbtion, daepfung, Cauchy1 cauchy2
                lb = [1 0]; ub = [20  100];
                options = optimoptions('fmincon');
                options = optimoptions(options,'Display', 'off');
                options = optimoptions(options,'FunValCheck', 'off');
                if tapecell
                    [x,fval,exitflag,output,lambda,grad,hessian] =...
                        fmincon(@retrieval_of_n_airref_thinlayer_tapecell,x,[],[],[],[],lb,ub,[],options);
                else
                    [x,fval,exitflag,output,lambda,grad,hessian] =...
                        fmincon(@retrieval_of_n_airref_thinlayer_pila,x,[],[],[],[],lb,ub,[],options);
                end
                global output_fit_opti
                S21DLN_targett(fitschleife) = output_fit_opti;
                nDLN(fitschleife)=x(1);
                kDLN(fitschleife)=x(2);
                S21DLN_processed(fitschleife) = S21target;
                mbranchpick(fitschleife)  = mmmmm;
                global testmich
                indexlauf(fitschleife)=testmich;
                %[nDLN(fitschleife),kDLN(fitschleife)] = retrieval_of_n_airref_thinlayer(nDLN(fitschleife-1),kDLN(fitschleife-1));
                nDLN_smooth(fitschleife) = (nDLN(fitschleife)+nDLN(fitschleife-1))/2;
                
            end
            nDLN_smooth(1)=nDLN_smooth(2);
            global npilaa
            global kpilaa
            npilaa(1) = 10; %please note that the first point is not calculated ...
            kpilaa(1) = 0;
            %% and more plotting
            if zwischenplott_fit
                figure(230+lauflauf)
                hold on
                plot(frequenzrange,abs(S21DLN_processed),'r')
                plot(frequenzrange,abs(S21DLN_targett),'k:')
                plot(frequenzrange,mbranchpick,'b')
                figure(260+lauflauf)
                hold on
                plot(frequenzrange,angle(S21DLN_processed),'r')
                plot(frequenzrange,angle(S21DLN_targett),'k:')
                figure(1021)
                hold on
                plot(frequenzrange,(nDLN),'-','color',farbelocal)
                plot(frequenzrange,(kDLN),':','color',farbelocal)
                figure(10021)
                plot(frequenzrange,indexlauf)
                
            end
        end
        %% and more calculating of refractive index
        if tapecell
            if fittenn
                npilaa = nDLN;
                kpilaa = kDLN;
            end
        else
            for fitschleife =2:length(frequenzrange)
                frequenzpunkt = frequenzrange(fitschleife);
                % refractive index of mixture and reference
                teflon =1;
                [npilaa(fitschleife),kpilaa(fitschleife)] = retrieval_of_n_airref_thinlayer_Mk2__Phasesep(npilaa(fitschleife-1),kpilaa(fitschleife-1));
                S21teflon_processed(fitschleife) = abs(S21target);
            end
            for fitschleife =2:length(frequenzrange)
                frequenzpunkt = frequenzrange(fitschleife);
                % refractive index of mixture and reference
                n_SiO2 = nteflon(fitschleife)+1i*kteflon(fitschleife);
                [sigma_ref_air(fitschleife),sigma_ref_SiO2(fitschleife)] = retrieval_of_sigma(n_SiO2);
            end
        end
        %% and plotting again, this time the conductivity...
        if sigmaplot
            if sigma_einzel
                figure(60+lauflauf)
                plot(frequenzrange,real(sigma_ref_air),'k')
                hold on
                plot(frequenzrange,imag(sigma_ref_air),'m')
                plot(frequenzrange,real(sigma_ref_SiO2),'b')
                plot(frequenzrange,imag(sigma_ref_SiO2),'c')
                legend('real air','imag air','real sio2','imag sio2')
            end
        end
        %% and more plotting
        if fittenn
            npilaa(1) = npilaa(2);
            kpilaa(1) = kpilaa(2);
            
            figure(105)
            hold on
            plot(frequenzrange,abs(S21DLN_processed),':k')
            if reference_on_air
                plot(freq,abs(S21prime),'r')
            else
                plot(freq,abs(S21DLN),'r')
            end
            xlim([0.5 3])
            title('comparing measured and calculated S21')
            ylabel('THz intens')
            xlabel('freq [THz]')
            figure(104)
            hold on
            plot(frequenzrange,unwrap(angle(S21DLN_processed)),':k')
            if reference_on_air
                plot(freq,unwrap(angle(S21prime)),'r')
            else
                plot(freq,unwrap(angle(S21DLN)),'r')
            end
            xlim([0.5 3])
            title('comparing measured and calculated S21')
            ylabel('phase')
            xlabel('freq [THz]')
            figure(106)
            plot(frequenzrange,(mbranchpick),'k')
            nDLN(1)=nDLN(2);
            kDLN(1)=kDLN(2);
        end
        %% exporting and substracting of reference samples if needed. If index2 > 2, this part is ignored
        temp_local = mod(temp,1000);
        if index2 ==1
            exportindex_F4 =exportindex_F4+1;
            exportindex = exportindex_F4;
            pilaf4_sigma_ref_SiO2(exportindex,1) = temp_local;
            pilaf4_sigma_ref_SiO2(exportindex,2:fitschleife+1) = sigma_ref_SiO2;
            pilaf4_sigma_ref_air(exportindex,1) = temp_local;
            pilaf4_sigma_ref_air(exportindex,2:fitschleife+1) = sigma_ref_air;
            if sigmaplot
                figure(1313)
                hold on
                plot(frequenzrange,(real(sigma_ref_air)),'color',farbelocal)
                plot(frequenzrange,(imag(sigma_ref_air)),':','color',farbelocal)
                sigi_loc = abs(real(sigma_ref_air))';
                freq_loc = frequenzrange';
            end
            save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_sigma.txt'],'sigi_loc','-ascii');
        end
        if index2==2
            exportindex_F3 =exportindex_F3+1;
            exportindex = exportindex_F3;
            pilaf3_sigma_ref_SiO2(exportindex,1) = temp_local;
            pilaf3_sigma_ref_SiO2(exportindex,2:fitschleife+1) = sigma_ref_SiO2;
            pilaf3_sigma_ref_air(exportindex,1) = temp_local;
            pilaf3_sigma_ref_air(exportindex,2:fitschleife+1) = sigma_ref_air;
            dif_sigma_ref_SiO2(exportindex,1)= temp_local;
            dif_sigma_ref_air(exportindex,1)= temp_local;
            %% calculates the conductivity from the permitivitty using the measured permitivitty of Aro mutant
            if lattice
                [argh indexxxx] = min(abs(pilaf4_sigma_ref_SiO2(:,1:2)-temp_local));
                dif_sigma_ref_SiO2(exportindex,2:fitschleife+1) = -sigma_ref_SiO2+pilaf4_sigma_ref_SiO2(max(indexxxx),2:end);
                dif_sigma_ref_air(exportindex,2:fitschleife+1) = -sigma_ref_air+pilaf4_sigma_ref_air(max(indexxxx),2:end);
            else
                dif_sigma_ref_air(exportindex,2:fitschleife+1) = sigma_ref_air;
                dif_sigma_ref_SiO2(exportindex,2:fitschleife+1) = sigma_ref_SiO2;
            end
            if sigmaplot
                figure(1008)
                hold on
                
                plot(frequenzrange,abs((real(dif_sigma_ref_SiO2(exportindex,2:fitschleife+1)))),':','color',farbelocal)
                plot(frequenzrange,abs((real(dif_sigma_ref_air(exportindex,2:fitschleife+1)))),'color',farbelocal)
                xlim([0.3 2])
            end
            sigi_loc = abs((real(dif_sigma_ref_air(exportindex,2:end))))';%abs((real(sigma__semianalytic)))';%
            freq_loc = frequenzrange';
            
            save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_sigma.txt'],'sigi_loc','-ascii');
            save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_frequenzrange_forallfits.txt'],'freq_loc','-ascii');
            if sigmaplot
                xlabel('frequency [THz]');
                ylabel('real(\sigma) [S/m]');
                title('real sigma air sigma_ref_SiO2  difference')
                legend(': is ref on SiO, surfaceimpedance',' is ref on air, surfaceimpedance','+numeric')
                figure(1009)
                hold on
                plot(frequenzrange,(imag(dif_sigma_ref_SiO2(exportindex,2:fitschleife+1))),':','color',farbelocal)
                plot(frequenzrange,(imag(dif_sigma_ref_air(exportindex,2:fitschleife+1))),'color',farbelocal)
                legend('ref on SiO, surfaceimpedance','ref on air, surfaceimpedance','numeric')
                xlabel('frequency [THz]');
                ylabel('imag(\sigma) [S/m]');
                title('imag sigma air sigma_ref_SiO2  difference')
            end
            figure(111)
            hold on
            if smooth
                plot(frequenzrange,real(nDLN_smooth),'color',farbelocal)
            else
                plot(frequenzrange,real(nDLN),'color',farbelocal)
            end
            %plot(frequenzrange,real(npilaa),'-.','color',farbelocal)
            if fittend
                % plot(frequenzrange,real(n_debye),':','color',farbelocal)
                dpilaaa(lauflauf) = dpila;
            end
            legend('real refractive, numeric','real refractive analytic')
            xlabel('frequency [THz]');
            ylabel('refractive index');
            title('refractive index of sample')
            figure(112)
            hold on
            plot(frequenzrange,imag(nDLN),'color',farbelocal)
            %plot(frequenzrange,imag(npilaa),'-.','color',farbelocal)
            if fittend
                %   plot(frequenzrange,imag(n_debye),':','color',farbelocal)
                legend('real refractive','debye model')
            end
            legend('real refractive numeric','imag refractive analytic')
            xlabel('frequency [THz]');
            ylabel('refractive index');
            title('refractive index of sample')
            
        end
        
        %% thats it kind of ...
        toc %time for one run is printed
        if fittenn %some more saving
            save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_npila.txt'],'nDLN','-ascii');
            save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_kpila.txt'],'kDLN','-ascii');
            save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_npila_numeric.txt'],'npilaa','-ascii');
            save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_kpila_numeric.txt'],'kpilaa','-ascii');
            save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_nteflon.txt'],'nteflon','-ascii');
            save([ordnerspeichern,probenname,'\',num2str(tempsaveing),'K_kteflon.txt'],'kteflon','-ascii');
        end
        if fittend
            dtotal(lauflauf) = dteflon
        end
    end
end
% the end




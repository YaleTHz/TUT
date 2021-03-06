
%% JN labbook 98,91,59-61,99,107-8,111,112,123-4, 2_20-4
%% General Remarks:
% All rights reserved (C) Jens Neu 2019, the software can be used according to GPL-3 licence.
% The software is provided as is without any warranty or reliability on my
% site. Please report bugs to jens.neu@yale.edu. Please do not distribute.
% uses on/off data, so you need a reference and the difference scan
% the first part defines some global variables, used to communicate with
% the fitting routines. The sample parameters need to be adjusted for
% relative volume and d_on for sample thickness or skin depth.
% If the sample is thicker than skin depth the n_substrat needs to be SiO2,
% if the sample is thinner I need to implement the etalon... needs to be done ... JN 12/08/17
% If the sample is sandwiched, the refractive index of air should be used
% as the refractive of the cover layer. 
% You also need a fairly accurate value for the non-photo excited
% refractive index. You can get that by measuring 2-5 accurately thick
% films on silicon. Silicon is better than SiO2, because the thickness has
% a smaller variation. 
% sno stands for the photo-layer, so if you look at nafion, it will still
% be called sno! 
% implemented Bruggeman to calculate the non excited layer parameters,
% using debye modell for Nafion
clear all %matlab says it is not needed .... usually ....
close all
tic %starts the internal glock, toc reads out
%% What you need to use this program:
% Reference measurement and difference measurement, (approx. line 103)
% Substrate measurements and non-excited measurements (loaded from  nmean_sno2 and nmean....
% functions: faxis.m, n_no_model_difference_v1.m,  n_DS_model_difference_v1.m import_for_timfiles.m importfiles_script.m
%% global variables
global Eon
global Eoff
global freq
global n_off
global n_air
global d_on
global d_penetration
global d_off
global d_gesamt
global Eon_durch_Eoff
global frequenzpunkt
global fschritt
global fmitteln
global n_substrat
global Etarget
global FractionE
global meff
global fmin
global fmax
global freq_sno
global f_sio2
global freq_long
global E_on_off_long
global fminDS
global fmaxDS
global x_fit
global y_fit
global sigmaa
global thinsample
global internal_FP
global fit_ds_transfer %exports fit function result from DS model for refractive index
global tapecell %parameter defining whether the sample is simply on substraet (0) or sandwiched between tape
global dteflon;
global cocker_DS
cocker_DS =0; %This parameter uses a different DS model.
%% Sample paramters
m_rel = 0.275; %electron Mass; Material Parameter only used for DS retrieval SnO2=0.275; Si =0.19;
d_on = 3.4775 *1e-6; % total film thickness also defined via layer numbers, in line 94
d_penetration = 10; % for 1*(0.105 *3/4)*1e-6; %value from [7]Green, M.A. and Keevers, M. "Optical properties of intrinsic silicon at 300 K ", Progress in Photovoltaics, p.189-92, vol.3, no.3; (1995)  the value is actually defined for bulk silicon, so not sure how well this applies to nanocores. the 3/4 adjustes for the volume fraction and is a  rought approx
nafion_literature =0; % Boolean, used to calculate the effective layer permitivitty from the literature values and Bruggeman, 
alpha =0; % fill factor for effective medium theory; from f_si:f_n = 2.8968:1
penetration_model =1; % boolean; in the case of a short penetration depth, this switch will replace the layer thickness with the penetration depth
internal_FP_wish =1; %  consideres the FP in the non-photo excited layer JN 2_23
tapecell =1; %if your sample is a sandwich sample, tape/material/tape
dteflon = 132/2; %thickness of one of the tape stripes
%% Physical constants, Only change in case of Maxwell's second coming
echarge = 1.602*1e-19; 
epsilon_0 = 8.854 *1e-12;
c0 =physconst('LightSpeed');
m_e = 9.1093829*1e-31; 
meff = m_e*m_rel; %effecitve mass for Drude Smith 
%% Load in data for the fourier transform (need time-domain files, .tim)
ordner= 'C:\ExperimentalData_backuped\brian\Brian\Conductive MOF\OPTP\'; %these two are building the load folder
measurementtype = ''; %this is a part of the folder name! 
ordnerspeichern = 'C:\ExperimentalData_backuped\brian\Brian\Conductive MOF\OPTP_results\'; %save folder
exportfolder =ordnerspeichern; %another varible for save folder ... 
exportname ='\samplename.txt'; % Part of the name of the generated file

token =4; %this is similar to the assigment sheet temprange picks 
if token ==1
onname = 'Dec26_001';
offname ='Dec26_002';
farbtoken ='r';
end
if token ==2
onname = 'Dec26_003';
offname ='Dec26_004';
farbtoken ='k';
end
if token ==3
onname = 'Dec26_005';
offname ='Dec26_006';
farbtoken ='b';
end
if token ==4
onname = 'Dec26_007';
offname ='Dec26_008';
farbtoken ='g';
end



layer =fix((token-2)/4)+1; % number of layers calculated from the position given in token. +1 accounts for counting from 1 instead of 0; /4 for the four lines per dataset
%% Excel export commands to allow the transport of the data to your collaborators, which do not use matlab and think that igor/origin is better ... Excel is fairly slow but the export function is need and more convinient than sending 100 .txt files to your collaborators
excelesheet = 'MOF'; % Really important that this is either Por, or RuP or Rup!! This varibale is later used to pick the correct penetration depth
excelesheetraw = [excelesheet ' Sample ' num2str(token)];
[time_on, amp_on] = importfiles_script([ordner,measurementtype,onname,'.tim']); %5L_D RuPhos  Sample
[time_off, amp_off] = importfiles_script([ordner,measurementtype,offname,'.tim']); %Ref
if tapecell
    if token ==1
        d_on = 146* 1E-6; % are these values correct....
    end
    if token ==2
        d_on = 163* 1E-6; % are these values correct....
    end
    if token ==3
        d_on = 138* 1E-6; % are these values correct....
    end    
    if token ==4
        d_on = 158* 1E-6; % are these values correct....
    end
    d_penetration = 15* 1E-6; 
else
    if layer ==1
        d_on =3.3* 1E-6;
    else
        d_on = (layer * 2.1 + 1.91) * 1E-6; % this value needs to be corrected as soon as we got the original data
    end
end
%% assigning penetration depth
if strcmp(excelesheet,'Por')
d_penetration = 1* 1E-6; 
end
if strcmp(excelesheet,'Rup')
d_penetration = 15* 1E-6;
end
if strcmp(excelesheet,'RuP')
d_penetration = 15* 1E-6;
end
%% dynamically picking the correct case depending on penetration length/sample thickness
% See "Applicability of the Thin-film Approximation in Terahertz Photoconductivity Measurements", J.Neu et al., APL 113 22 2018
if d_on < d_penetration %defines whether the thin layer assumption is used and the substrat is the substrat, or if the layer is seen as the full layer
    thinsample =0; 
    internal_FP =0; %when the full layer is photo-excited, there can't be a FP in the non-excited, because there is no non-excited
    d_thinfilmapprox = d_on;
    d_gesamt = d_on; 
else
    thinsample =1;
    d_thinfilmapprox = d_penetration;
    internal_FP =0; %without the penetration model, the full layer is assumed conducting, so not non-excited layer, no FP in this layer
    if penetration_model
        d_gesamt = d_on; %historically the layer was assumed to be totally excited... this redefinition safes the thickness in a new variable
        d_off = d_on - d_penetration;
        d_on = d_penetration; %not sure about that .... it will reduce the layer thickness to the values for the penetration depth. So assuming that all charges are in this range. This is needed to calculate quatum efficiency
        internal_FP =internal_FP_wish; % The internal FP are considered if you want them to be considered... 
    end
end

%% Settings for fitting and plotting 
fmitteln=1; %averaging the values in between the steps
DS =1; %DS fitting with and without a model 
zwischenplotten =1;% intermideate ploits 
thinfilm =1; %using also thin film approx. 
DSfitting_after_retrieval = 1; %fitting DS after the fact on the data, not yet implemented, needs to be on for export
fmin =0.25; fmax =3;fschritt =0.01; %make sure that the loaded un-excited data has a sufficient spacing and size.
fminDS = 0.5; fmaxDS =2.2; % Freq range of DS model to calculatre n
fmin_DS_fit =0.3; fmax_DS_fit =2.5; % for the fitting range after calculation, the cases enforce that this is smaller than the full range. 
bruggeman =0; % Boolean, calculates back from the layer refractive index to the inclusion, using Bruggeman-equations
exportieren =1; %defines whether it should write into an excel file or not. 
speicherort =ordnerspeichern; %defines folder of excel file 
%% some stuff

f_N = 1/(1+alpha); f_Si = 1 - f_N; %volum fractions for Bruggeman
sigi = '\sigma [S/m]'; %defines labels for the figures
freqlabel ='frequency [THz]';%defines labels for the figures
MS =8; %defines markersize for plotting
%% Take the fourier transform.
NNT = 2^18; %zero padding
freq_on = faxis(time_on,NNT)';
freq_off = faxis(time_off,NNT)';
fft_on = conj(fft(amp_on,NNT))'; %difference sno2 photo
fft_off = conj(fft(amp_off,NNT))'; %reference of sno2 non photo
freqdiff = freq_on-freq_off;
testvalue = sum(abs(freqdiff)); %if test value is larger than zero you know that you messed soemthign up 
freq= freq_on;
freq_long =freq_on;

Eon_test = (fft_on + fft_off); %Building the complex Eon value; the postfix _test is used for all original data, the data used for calculation is reduced to 561 points, because that is the resultion of the n_substrat and n_off; need to be adjusted if different data point resolution is used
Eoff_test = fft_off;
E_on_off_long = Eon_test./Eoff_test; % long is the same as test
%% input from previous meassurements
if tapecell %only works if both refractive index have same spacing, so far ... 
    main_previous_results = 'C:\ExperimentalData_backuped\brian\Brian\Conductive MOF\TDS_results\';
    nmean_quartz = load([main_previous_results,'mean_nmean']); %loads the frequency dep. index of tape
    kmean_quartz = load([main_previous_results,'mean_kmean']); %loads the frequency dep. imag. index of tape
    nmean_sno2_test = load([main_previous_results,'mean_n_MOF_mean']); %loads parameters for MOF
    kmean_sno2_test = load([main_previous_results,'mean_k_MOF_mean']); %loads parameters for MOF
    freq_sno_test = load([main_previous_results,'freq']); %load frequency of MOF
    f_sio2 = load([main_previous_results,'freq']);%freq of sio2
else
    nmean_quartz = load('nmean.txt'); %loads the frequency dep. index of quartz
    kmean_quartz = load('kmean.txt'); %loads the frequency dep. imag. index of quartz
    nmean_sno2_test = load('nmean_sno2'); %loads parameters for SnO2 on silicon
    kmean_sno2_test = load('kmean_sno2'); %loads parameters for SnO2 on silicon
    freq_sno_test = (load('freq_sno.txt')); %load frequency of sno2
    f_sio2 = (load('freq.txt'));%freq of sio2
end
freq_size_test = size(f_sio2) - size(freq_sno_test);
n_size_test = size(nmean_sno2_test) - size(nmean_quartz);
argh = [freq_size_test n_size_test];

%% calculation the refractive index from literature, alternative to the previous input case
if nafion_literature %[5]Yurchenkoa and Zaytsevb, JAP, 116, 113508 (2014) https://doi.org/10.1063/1.4896194
    eps_inf = 2.2; delta_eps1 = 1.17; tau1= 7.2; delta_eps2 = 0.44; tau2 = 0.12;
    eps_lit = eps_inf + (delta_eps1./(1+ 2.*pi.*1i.*freq_sno_test.*tau1))+ (delta_eps2./(1+ 2.*pi.*1i.*freq_sno_test.*tau2)); %in THz units, NOT SI, be careful when editing
    % the imaginary part has a different sign than in [5]!
    n_Si = 3.4178.*((freq_sno_test*0)+1);% [6] Grischkowsky, JOSA B, 7,10 (1990) Grishttp://josab.osa.org/abstract.cfm?URI=josab-7-10-2006
    eps_Si = n_Si.^2; % The previous formula assumes loss-free silicon, which is a decent approximation
    eps_lit = conj(eps_lit); % accounting for the previous sign error in the debye model they used.... but in the paper the simply flipped it and ignored that they are using the electrical engineering formula....
    %% calculating the nafion film via Bruggeman
    %eq 29 from [7] A. Markel, JOSA A, 33,7 (2016) http://dx.doi.org/10.1364/JOSAA.33.001244
    b = (2*f_Si - f_N) .* eps_Si +((2*f_N - f_Si) .* eps_lit);
    eps_layer = (b + sqrt(8*eps_Si.*eps_lit + b.^2))./4;    
    n_comp = sqrt(eps_layer);
    nmean_sno2_test = real(n_comp);
    kmean_sno2_test = imag(n_comp); %accounting for the previous sign error in the debye model they used.... but in the paper the simply flipped it and ignored that they are using the electrical engineering formula....
end
%% ensuring that refractifve index have the same length, only works if they
%have already the same spacing! 
if max(freq_sno_test) > max(f_sio2)
    maxfreq_ref = find(freq_sno_test>max(f_sio2),1,'first');
    freq_sno = freq_sno_test(1:maxfreq_ref-1);
    nmean_sno2 = nmean_sno2_test(1:maxfreq_ref-1);
    kmean_sno2 = kmean_sno2_test(1:maxfreq_ref-1);
else
    freq_sno = freq_sno_test;
    nmean_sno2 = nmean_sno2_test;
    kmean_sno2 = kmean_sno2_test;
end
%Interpolate lines for the index of refraction
n_off = (nmean_sno2 - 1i*kmean_sno2);
size(n_off)

if thinsample %wether the substrate is the non photoexcited or the fused silica
    n_substrat = n_off;
else
    n_substrat = (nmean_quartz - 1i*kmean_quartz);
end
n_sio2 = (nmean_quartz - 1i*kmean_quartz);
maxfreq_exp = find(freq>(max(f_sio2)),1,'first');
Eon = Eon_test;
Eoff = Eoff_test;
%Using S parameter for Simulation, For experiment we can use the experimental Data for the two on/off measurements Labbook JN 123-124
clear freq
%% Discretization of input data, limit it to 0.2 to 3 THz, with a 0.005 sampling
for lauf = 1:length(freq_sno)
    freqeuncy_loc = freq_sno(lauf); 
    if fmitteln
        index_lb = find(freqeuncy_loc>(freq_long+fschritt/2),1,'last');
        index_ub = find(freqeuncy_loc>(freq_long-fschritt/2),1,'last');
        Eon_durch_Eoff(lauf) = mean(E_on_off_long(index_lb:index_ub));
        Eon_short(lauf) = mean(Eon(index_lb:index_ub));
        Eoff_short(lauf) = mean(Eoff(index_lb:index_ub));        
        freq(lauf) = mean(freq_long(index_lb:index_ub));
    else
    end
end
testtest = sum(abs(freq-freq_sno)); % if this value is larger 0, something is wrong; If this value is smaller than 0, something is terribly wrong
%% checking if the discretized length vectors and the original "look" the same
n_air = 1+(n_substrat*0); % ensures that the refractive index vector of air is the same length and freq range as SiO2
% indexstart = find(freq>fmin,1,'first');
% indexend = find(freq>fmax,1,'first');
% freq = freq(1,indexstart:indexend);
% Eref = Eref(1,indexstart:indexend);
% Eon = Eon(1,indexstart:indexend);
if zwischenplotten % 
   figure(5)
   plot(freq_long,abs(E_on_off_long),'k');
   title('difference signal, full and truncated')
   hold on
   plot(freq,abs(Eon_durch_Eoff),'r+');
   xlim([fmin fmax+1])
   figure(6)
   plot(freq_long,abs(Eon),'k');
   hold on
   plot(freq,abs(Eon_short),'r+');
   title('Eon signal, full and truncated')
   xlim([fmin fmax+1])
   figure(7)
   plot(freq_long,abs(Eoff),'k');
   hold on
   plot(freq,abs(Eoff_short),'r+');
   xlim([fmin fmax+1])
   title('Eoff signal, full and truncated')
   figure(8)
   plot(time_on,amp_on)
end
%% Retrieval of n, no model
eps_nafion(1)=1;eps_Si_reduced(1)=1; % initializing
frequenzrange = fmin:fschritt:fmax;
n_on(1) = n_off(1); %initialization of fit 
lb = [1 0]; ub = [1000 1000]; %boundaries
for fitschleife =2:length(frequenzrange) 
    frequenzpunkt = frequenzrange(fitschleife);
    x0(1) = real(n_on(fitschleife-1)); %intialization of the n-th fit point alternatively can also intialize with the n_off values; as long as the change is small, would be fine
    x0(2) = imag(n_on(fitschleife-1));
    options = optimoptions('fminimax');
    options = optimoptions(options,'Display', 'off');
    [x,~,~,exit_nomodel_opti,~,~] = ...
        fminimax(@new__Eon_durch_Eoff_no_model_v0,x0,[],[],[],[],lb,ub,[],options); %See Labbook 123 JN, tape cell implemented
    n_on(fitschleife) = x(1) - 1i*x(2); %export out of loop
    if thinsample % Thin layer have rather large refractive index, so if the initializing of the algorithm is not good, the boundaries need to be large to allow convergence
    else
        lb = [x(1)*0.6 x(2)*0.6]; % avoiding 2pi jumps and unphysical fast changes.
        ub = [x(1)*1.8 x(2)*1.8];
    end
    Etarget_compl(fitschleife) =   Etarget; %For plotting
    FractionE_compl(fitschleife) =    FractionE; %For plotting
    %% reduces the size of the input epsilons to the output size, this is needed for EMT calculations at line 390
    lower_limit = find(freq_sno_test>(frequenzpunkt-(fschritt/2)),1,'first');
    upper_limit =find(freq_sno_test>(frequenzpunkt+(fschritt/2)),1,'first');
    if bruggeman
        eps_nafion(fitschleife) = mean(eps_lit(lower_limit:upper_limit));
        eps_Si_reduced(fitschleife) = mean(eps_Si(lower_limit:upper_limit));
    end
end
n_on(1) = n_on(2);
if bruggeman
    eps_nafion(1)=eps_nafion(2);
    eps_Si_reduced(1) = eps_Si_reduced(2);
end
%if zwischenplotten %plotting the goal and the calculated function of the fit
    figure(42)
    hold on
    plot(frequenzrange,angle(FractionE_compl),'r','LineWidth',2)
    plot(frequenzrange,angle(Etarget_compl),'k:','LineWidth',3)
    title('phase optimization, mess and calc')
    figure(21)
    hold on
    plot(frequenzrange,abs(FractionE_compl),'r','LineWidth',2)
    plot(frequenzrange,abs(Etarget_compl),'k:','LineWidth',3)
    title('abs optimization, mess and calc')
%end
clear x
clear lb
clear ub
%% Retrieval of n, Drude Smith
if DS
    %% Drude parameter, See 10 lines below for units and what is what.
    x0 =[0.01 30 -0.5];  
    lb =[0.00001 0.1 -1];
    ub =[inf inf 0];
    options = optimoptions('fminimax');
    %% Modify options setting
    options = optimoptions(options,'Display', 'off','FiniteDifferenceType','central');
    [x,~,~,exit_DS_opti,output_otpimization_stuff,~] = ...
        fminimax(@new__Eon_durch_Eoff_DS_v0,x0,[],[],[],[],lb,ub,[],options); %supports tapecell
    omega = 2*pi*freq*1e12;
    k0 = omega/c0;
    N = x(1)*1e23; %in 10^17/(cm*cm*cm)
    tau = x(2)*1e-15; %in fs
    c = x(3); %c parameter
    x_from_model =x;
    %% calculating for plotting
    epsilon_l = n_off.^2; %lattice parameter; maybe times 0+1 ...
    dn = 1- 1i*omega*tau;
    omegapp = (N.*echarge.*echarge.*tau./meff);
    sigma = omegapp./(dn) .*(1+ c./(dn));
    epsilon_on = epsilon_l + ((1i*sigma)./(epsilon_0*omega));
    n_on_DS = sqrt(conj(epsilon_on));
    DS_DS_model=x;
    figure(13)
    plot(freq,abs(Eon_durch_Eoff),'k')
    hold on
    plot(freq,abs(fit_ds_transfer),'r')
    figure(14)
    plot(freq,angle(Eon_durch_Eoff),'k')
    hold on
    plot(freq,angle(fit_ds_transfer),'r')

end
renorm = (fschritt/0.005*10);
n_off_scale = n_off(1:renorm:end);
f_off_scale = freq_sno(1:renorm:end); 
%% Calculate conudctivity of the sample, from the refractive index, without model

conductivity = -1.*1i.*f_off_scale*2*pi*1e12.*((n_on.^2).*epsilon_0 - (n_off_scale.^2).*epsilon_0);
%% Comming soon: Fitting DS on sigma which has been calculated without model
if DSfitting_after_retrieval
    if fmin > fmin_DS_fit
        fmin_DS_fit = fmin;
    end
    if fmax < fmax_DS_fit
        fmax_DS_fit = fmax;
    end
    anfangindex=find(fmin_DS_fit<f_off_scale,1,'first')-1;
    endindex=find(fmax_DS_fit<f_off_scale,1,'first')-1;
    frequenzrange_fit = fmin_DS_fit:fschritt:fmax_DS_fit;
    x_fit = frequenzrange_fit;
    y_fit = -conj(conductivity(anfangindex:endindex));
    XXX = frequenzrange_fit;
    YYY = -conj(conductivity(anfangindex:endindex));
    
    % Drude parameter
    if DS
        x0 = [x(1) x(2) x(3)];
    else
        x0 = [1 60 0];
    end
    
    lb =[0.00001 0.1 -1];
    ub =[inf inf 0];
    options = optimoptions('fminimax');
    % Fitting
    options = optimoptions(options,'Display', 'off');
    [x,fval,maxfval,exitflag,output,lambda] = ...
        fminimax(@DSfit,x0,[],[],[],[],lb,ub,[],options);
    
    DS_ergebnis_optitool= x ;
    %% lsq fit
    v =x_from_model;
    ydata = [real(YYY)',imag(YYY)'];
    %objfcn = @(x,xdata)1.0245*1e-9 .*((x(1)*x(2))./(1-1i*xdata*2*pi*1e12*x(2)*1e-15) .*(1+ (x(3)./(1-1i*xdata*2*pi*1e12*x(2)*1e-15)))          );
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
    % DS_seperate
    lb =[0.0001 0.5 -1];
    ub =[inf inf 0];
    %   lb =[]; ub=lb;
    %  x = lsqcurvefit(DS_seperate,v,XXX,YYY,lb,ub,options);
    [vestimated,resnorm,resid,exitflag,output,lambda,J] = ...
        lsqcurvefit(@DS_seperate,v,XXX,ydata,lb,ub,options);
    dsfitresult=DS_seperate(vestimated,frequenzrange);
    dsfit_complex = dsfitresult(:,1)+1i*dsfitresult(:,2);
    %DSfit(x_from_model);
    DS_ergebnis_fit= vestimated;
    DS_ergebnis_ncalc= x_from_model;
    ci = nlparci(vestimated,resid,'jacobian',J);
    DS_error_fit = abs(ci(:,1)-ci(:,2));
    oneterm = sum(sum((ydata-mean(ydata)).^2));
    R2 = 1- (resnorm/oneterm);
    DS_resnorm = R2;
end
%% export need to be defined better 
zexport = [(f_off_scale(8:37).*1E12)' -1*real(conductivity(8:37))' imag(conductivity(8:37))'];
save([exportfolder,exportname],'zexport','-ascii')
%% Thin Film approximation, 
if thinfilm
    Z0 = 377;
    sindex = real(n_substrat);
    W_txMeas = abs(Eon_short)./abs(Eoff_short);
    diff = angle(Eon_short);
    ref = angle(Eoff_short);
    diff_un = unwrap(diff);
    ref_un = unwrap(ref);
    W_argMeas_angle = (diff_un-ref_un);
    W_argMeas = atan2(-imag(Eon_short./Eoff_short),real(Eon_short./Eoff_short));
    W_argMeas_un = unwrap(W_argMeas);
    cons = (sindex + 1)./(Z0.*d_thinfilmapprox);
    W_realCond = (cons).*((cos(W_argMeas)./W_txMeas)- 1);
    W_imagCond =-1*(cons).*(sin(W_argMeas)./W_txMeas);
    W_realCond_angle = (cons).*((cos(W_argMeas_angle)./W_txMeas)- 1);
    W_imagCond_angle =(cons).*(sin(W_argMeas_angle)./W_txMeas);
    trans = Eon_short./Eoff_short;
    cond = cons.*((1./trans)-1);
    r_cond = real(cond);
    i_cond = imag(cond);
    sigma_thinfilm = r_cond - 1i*i_cond;
    if DSfitting_after_retrieval
        v =x_from_model;
        ydata = [real(sigma_thinfilm)',imag(sigma_thinfilm)'];
        %objfcn = @(x,xdata)1.0245*1e-9 .*((x(1)*x(2))./(1-1i*xdata*2*pi*1e12*x(2)*1e-15) .*(1+ (x(3)./(1-1i*xdata*2*pi*1e12*x(2)*1e-15)))          );
        options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
        % DS_seperate
        lb =[0.0001 0.5 -1];
        ub =[inf inf 0];
        %   lb =[]; ub=lb;
        %  x = lsqcurvefit(DS_seperate,v,XXX,YYY,lb,ub,options);
        [vestimated,resnorm,resid,exitflag,output,lambda,J] = ...
            lsqcurvefit(@DS_seperate,v,freq,ydata,lb,ub);
        dsfitresult_thin=DS_seperate(vestimated,freq);
        dsfit_thinfilm = dsfitresult_thin(:,1)+1i*dsfitresult_thin(:,2);
        DS_ergebnis_fit_thinfilm= vestimated;
        ci = nlparci(vestimated,resid,'jacobian',J);
        DS_error_fit_thinfilm = abs(ci(:,1)-ci(:,2));
                oneterm = sum(sum((ydata-mean(ydata)).^2));
        R2 = 1- (resnorm/oneterm);
        resnorm_thinfilm = R2;
    end
end
%% Plotting of result
figure(26) 
plot(f_off_scale, -1*real(conductivity),farbtoken)
hold on
plot(f_off_scale, imag(conductivity),farbtoken,'LineStyle',':')
figure(24) 
plot(f_off_scale, -1*real(conductivity),'k')
hold on
plot(f_off_scale, imag(conductivity),'k:')
lelele{1} = 'real \sigma, NM';
lelele{2} = 'imag \sigma, NM';
if DS
    plot(freq(1:10:end),real(sigma(1:10:end)),'r')
    plot(freq(1:10:end), imag(sigma(1:10:end)),'r:')
    lelele{end+1} = 'real \sigma, DS as model';
    lelele{end+1} = 'imag \sigma, DS as model';
end

if thinfilm
    plot(freq,real(sigma_thinfilm),'b')
    plot(freq,imag(sigma_thinfilm),'b:')
    lelele{end+1} = 'real \sigma, thinfilm';
    lelele{end+1} = 'imag \sigma, thinfilm';
    figure(27)
        plot(freq,real(sigma_thinfilm),farbtoken)
        hold on
        plot(freq,imag(sigma_thinfilm),farbtoken,'LineStyle',':')
        figure(24)
    if DSfitting_after_retrieval
        plot(freq,real(dsfit_thinfilm),'c')
        plot(freq,imag(dsfit_thinfilm),'c:')
        lelele{end+1} = 'real \sigma, thinfilm DSfit';
        lelele{end+1} = 'imag \sigma, thinfilm DSfit';     
    end
end

if DSfitting_after_retrieval %comming soon
    plot(frequenzrange,real(dsfit_complex),'m')
    plot(frequenzrange,imag(dsfit_complex),'m:')
    plot(x_fit,real(sigmaa),'g')% from DS fit after model
    plot(x_fit,imag(sigmaa),'g:')   
    lelele{end+1} = 'real \sigma, DS fit after';
    lelele{end+1} = 'imag \sigma, DS fit after';
    lelele{end+1} = 'real \sigma, DS fit after';
    lelele{end+1} = 'imag \sigma, DS fit after';
end





figure(24)
legend(lelele,'Location','EastOutside')
title('Conductivity of the film')
xlabel(freqlabel)
ylabel(sigi)

%igor_thinfilm = [freq' real(sigma_thinfilm)' imag(sigma_thinfilm)'];
%igor_extensive = [f_off_scale' -1*real(conductivity)' imag(conductivity)'];


figure(1)
plot(frequenzrange(2:end),real(n_on(2:end)),'*r','MarkerSize',MS) %(2:end) because the first point is the intialization value
hold on
plot(frequenzrange(2:end),-imag(n_on(2:end)),'*b','MarkerSize',MS)
plot(freq_sno(2:end),real(n_off(2:end)),'r','MarkerSize',MS) %(2:end) because the first point is the intialization value
hold on
plot(freq_sno(2:end),-imag(n_off(2:end)),'b','MarkerSize',MS)
title('refractive index, no model and DS model for retrieval')
if DS
    plot(freq(1:10:end),real(n_on_DS(1:10:end)),'m','linewidth',1.5)
    plot(freq(1:10:end),-imag(n_on_DS(1:10:end)),'c','linewidth',1.5)
end

%igordata = [ 

legend('real(no model on)','imag(no model on)','real(no model off)','imag(no model off)','real(DS model on)','imag(DS model on)')
%plot(freqeps,real(n_target),'g:','linewidth',2.0)
%plot(freqeps,-imag(n_target),'c:','linewidth',2.0)
%legend('n on, no model, real','n on, no model, imag','n on, DS, real','n on, DS, imag','n on input, real','n on input, imag')
xlabel(freqlabel)
ylabel('refractive index')
%% Approximation of Etalon in the sample, is calculated only for thin samples. JN p128
% this value can be used after the fact to justify ignoring the etalon in
% the non photo excited sample. 
if thinsample
    if internal_FP_wish
    else
        rho_on_off = (n_off(1:10:end) - n_on)./(n_off(1:10:end) + n_on);
        rho_off_sub = (n_off(1:10:end) - n_sio2(1:10:end))./(n_off(1:10:end) + n_sio2(1:10:end));
        figure(100)
        hold on
        plot(frequenzrange,abs(rho_on_off),'k')
        plot(frequenzrange,abs(rho_off_sub),'r')
        plot(frequenzrange,abs(rho_off_sub.*rho_on_off),'b')
        testvalue_internalref = max(abs(rho_off_sub.*rho_on_off))*100
        legend('reflection on off','reflection off substrate','reflection both')
        title('chekcs error from ignoring internal FP')
        figure(101)
        plot(frequenzrange,abs(rho_off_sub.*rho_on_off)*100,'b')
        title('chekcs error from ignoring internal FP')
        if testvalue_internalref > 1
            dim = [.2 .5 .3 .3];
            str = ['ARGH! the effect of the etalon in the nonswitchted substrat is', newline,', getting larger  ... would be/is ',num2str(testvalue_internalref),' %'];
            annotation('textbox',dim,'String',str,'FitBoxToText','on');
        end
    end
end
if bruggeman
    eps_on = n_on.^2; %value from numerical calculation of the layer refractive index
    nominator = 2.*eps_on.*eps_on + (eps_nafion.*eps_on.*(f_Si - 2*f_N));
    denominator = eps_nafion + eps_on.*(2*f_Si - f_N);
    eps_si_on = conj((nominator)./(denominator));
    delta_eps = eps_si_on-eps_Si_reduced;
    omegaa = frequenzrange *1E12*2*pi;
    sigma_no_model_EMT = -1i.*omegaa.*delta_eps.*epsilon_0;
    if thinfilm
        sigma_thinfilm_EMT = -(sigma_thinfilm./(0.5*f_N-f_Si));
    end
    
    figure(1001)
    hold on
    title('compares sigma, after EMT')
    plot(frequenzrange,real(sigma_no_model_EMT),'k')
    plot(frequenzrange,imag(sigma_no_model_EMT),'k')
    lalelu{1} = 'real sigma EMT no model';
    lalelu{2} = 'imag sigma EMT no model';
    if DSfitting_after_retrieval
        y_fit = (sigma_no_model_EMT);
        x_fit = frequenzrange;
        x0 =DS_DS_model;
        %% Start with the default options
        options = optimset;
        %% Modify options setting
        options = optimset(options,'Display', 'off');
        [x,fval,exitflag,output] = ...
            fminsearch(@DSfit,x0,options);
        plot(frequenzrange,real(sigmaa),'r')
        plot(frequenzrange,imag(sigmaa),'r+')
        lalelu{end+1} = 'real sigma DS no model';
        lalelu{end+1} = 'imag sigma DS no model';
        DS_EMT_nomodel=x;
    end
    if thinfilm
        plot(f_sio2,real(sigma_thinfilm_EMT),'b')
        plot(f_sio2,imag(sigma_thinfilm_EMT),'b-.')
        lalelu{end+1} = 'real sigma thinfilm approx.';
        lalelu{end+1} = 'imag sigma thinfilm approx.';
        if DSfitting_after_retrieval
        y_fit = (sigma_thinfilm_EMT);
        x_fit = f_sio2;
        x0 =DS_DS_model;
        %% Start with the default options
        options = optimset;
        %% Modify options setting
        options = optimset(options,'Display', 'off');
        [x,fval,exitflag,output] = ...
                fminsearch(@DSfit,x0,options);
            plot(f_sio2,real(sigmaa),'c')
            plot(f_sio2,imag(sigmaa),'c-.')
            lalelu{end+1} = 'real sigma DS thinfilm';
            lalelu{end+1} = 'imag sigma DS thinfilm';  
            DS_EMT_thinfilm =x;
        end
    end


    legend(lalelu,'Location','EastOutside')
    box on
    %% not yet working correct
    xlabel(freqlabel)
    ylabel(sigi)
end

%d_thinfilmapprox % thickness in which the chrages are
Nphoton_total = 4.1489e+19; %per m^2  calculated for knifeedge measurement at 150mW pump power,C:\ExperimentalData_backuped\Jens\thinfilm\knifeedgeplotter.m
%% calculating the total conversion efficiency per m^2
if DS
N_e_density_DS_ncalc = DS_ergebnis_ncalc(1);
N_e_DS_ncalc_per_m2 = N_e_density_DS_ncalc*d_thinfilmapprox*1e23; %the 23 is because I pulled that value out of the fit
qeff_total_DS_ncalc = N_e_DS_ncalc_per_m2/Nphoton_total;
end
if DSfitting_after_retrieval
N_e_density_DS = DS_ergebnis_fit(1); 
N_e_DS_per_m2 = N_e_density_DS*d_thinfilmapprox*1e23;
qeff_total_DS = N_e_DS_per_m2/Nphoton_total;

N_e_density_DS_thin = DS_ergebnis_fit_thinfilm(1);
N_e_DS_thin_per_m2 = N_e_density_DS_thin*d_thinfilmapprox*1e23;
qeff_total_DS_thin = N_e_DS_thin_per_m2/Nphoton_total;

N_e_density_DS_optitool = DS_ergebnis_optitool(1);
N_e_DS_optitool_per_m2 = N_e_density_DS_optitool*d_thinfilmapprox*1e23;
qeff_total_DS_optitool = N_e_DS_optitool_per_m2/Nphoton_total;
end
if exportieren & DSfitting_after_retrieval & DS
    if tapecell
                figure(24)
        drudeparameter =[DS_ergebnis_ncalc DS_ergebnis_fit DS_ergebnis_fit_thinfilm DS_ergebnis_optitool 0 0 ...
            DS_error_fit' DS_error_fit_thinfilm' 0 0 ...
            N_e_DS_ncalc_per_m2/1e17 N_e_DS_per_m2/1e17 N_e_DS_thin_per_m2/1e17 N_e_DS_optitool_per_m2/1e17 0 0 ...
            qeff_total_DS_ncalc*100 qeff_total_DS*100 qeff_total_DS_thin*100 qeff_total_DS_optitool*100 0 0 ...
            DS_resnorm/10000 resnorm_thinfilm/10000];
        %speicherort ='C:\ExperimentalData_backuped\Jens\thinfilm\evaluation_and_fit\';
        save([speicherort,onname,'.mat']);
        savefig([speicherort,onname,'.fig']);
        buchstabe ='B';
        xlswrite([speicherort,'Overview.xlsx'],drudeparameter,[excelesheet, ' fit results'],[buchstabe,num2str(token+2)]);
        xlswrite([speicherort,'Overview.xlsx'],onname,[excelesheet, ' fit results'],['AM',num2str(token+2)]); % make sure it actually points to an empty cell
        
        sigicompl = [f_off_scale' real(conductivity)' imag(conductivity)'];
        xlswrite([speicherort,'Overview.xlsx'],sigicompl,excelesheetraw,'A2');
        if thinfilm
            sigicompl = [freq(1:10:end)' real(sigma_thinfilm(1:10:end))' imag(sigma_thinfilm(1:10:end))'];
            xlswrite([speicherort,'Overview.xlsx'],sigicompl,excelesheetraw,'E2');
        end
        if DS
            sigicompl = [freq(1:10:end)' real(sigma(1:10:end))' imag(sigma(1:10:end))'];
            xlswrite([speicherort,'Overview.xlsx'],sigicompl,excelesheetraw,'I2');
        end
        
    else
        figure(24)
        drudeparameter =[DS_ergebnis_ncalc DS_ergebnis_fit DS_ergebnis_fit_thinfilm DS_ergebnis_optitool 0 0 ...
            DS_error_fit' DS_error_fit_thinfilm' 0 0 ...
            N_e_DS_ncalc_per_m2/1e17 N_e_DS_per_m2/1e17 N_e_DS_thin_per_m2/1e17 N_e_DS_optitool_per_m2/1e17 0 0 ...
            qeff_total_DS_ncalc*100 qeff_total_DS*100 qeff_total_DS_thin*100 qeff_total_DS_optitool*100 0 0 ...
            DS_resnorm/10000 resnorm_thinfilm/10000];
        %speicherort ='C:\ExperimentalData_backuped\Jens\thinfilm\evaluation_and_fit\';
        save([speicherort,onname,'.mat']);
        savefig([speicherort,onname,'.fig']);
        buchstabe ='G';
        xlswrite([speicherort,'Overview.xlsx'],drudeparameter,excelesheet,[buchstabe,num2str(token)]);
        xlswrite([speicherort,'Overview.xlsx'],onname,excelesheet,['AQ',num2str(token)]); % make sure it actually points to an empty cell
    end
end




toc


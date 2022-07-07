close all
clear variables

flags.paramset = 100; %choose the parameters to run the model
%paramset
% -1 - artificial parameters for non-dimensional model
% 0 - population parameters from Chapin model
% 1 - Pfizer parameters
% 2 - Moderna parameters
% 100 - these are the different data sets from the literature given by

flags.ID = 3;
flags.save = 1; %to save (1) or not (0)
%% 

run_params;

scales = [1;params.V0;params.T0;params.B0;params.A0;params.C0;params.F0;params.I0];

t_final = 365;
dose_times = linspace(1,t_final,t_final);

y0 = [1;0;0;0;0;0;0;0];

odeopts = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-9);

flags.modtype='scaled';

params.alphaL = 0;
params.lambdaC=0;
params.lambdaF=0;
% params.kappaI=0;

peak = zeros(length(dose_times)+1,1);
peak_time = zeros(length(dose_times)+1,1);
tot_imm = zeros(length(dose_times)+1,1);
FWHM_spot1 = zeros(length(dose_times)+1,1);
FWHM_spot2 = zeros(length(dose_times)+1,1);
delT = zeros(length(dose_times)+1,1);

[t,y] = ode23t(@ODEf,[0,t_final/params.t0],y0./scales,odeopts,params,flags);

peakspot = find(y(:,5)==max(y(:,5)));
peak(1) = y(peakspot,5);
peak_time(1) = t(peakspot);
tot_imm(1) = trapz(t,y(:,5));
FWHM_spot1(1) = find(y(:,5)>peak(1)/2,1);
FWHM_spot2(1) = FWHM_spot1(1)-1 + find(y(FWHM_spot1(1):end,5)<peak(1)/2,1);
delT(1) = t(FWHM_spot2(1))-t(FWHM_spot1(1));
A{1} = y(:,5);
Time{1} = t;

parfor k=1:length(dose_times)
    t_dose = dose_times(k);
    sol = ode23t(@ODEf,[0,t_dose/params.t0],y0./scales,odeopts,params,flags);
    y02 = sol.y(:,end);
    y02(1) = y02(1)+1;
    sol = odextend(sol,[],t_final/params.t0,y02,odeopts,params,flags);
    t = sol.x';
    y = sol.y';
    
%     Ls = y(:,1);
%     Vs = y(:,2);
%     Ts = y(:,3);
%     Bs = y(:,4);
%     As = y(:,5);
%     Cs = y(:,6);
%     Fs = y(:,7);
%     Is = y(:,8);
    peakspot = find(y(:,5)==max(y(:,5)));
    peak(k+1) = y(peakspot,5);
    peak_time(k+1) = t(peakspot);
    tot_imm(k+1) = trapz(t,y(:,5));
    FWHM_spot1(k+1) = find(y(:,5)>peak(k+1)/2,1);
    FWHM_spot2(k+1) = FWHM_spot1(k+1)-1 + find(y(FWHM_spot1(k+1):end,5)<peak(k+1)/2,1);
    delT(k) = t(FWHM_spot2(k+1))-t(FWHM_spot1(k+1));
    A{k+1} = y(:,5);
    Time{k+1} = t;
end

figure
plot(dose_times,peak(2:end)*params.A0,'linewidth',2);
figure
plot(dose_times,peak_time(2:end)*params.t0,'linewidth',2);
figure
plot(dose_times,peak_time(2:end)-dose_times','linewidth',2)
figure
plot(dose_times,tot_imm(2:end)*params.t0*params.A0,'linewidth',2)
figure
plot(dose_times,delT(2:end)*params.t0,'linewidth',2);

for k=1:length(dose_times)
    temp = sign(diff(A{k}));
    spots = find(temp==0);
    if ~isempty(spots)
        for j=1:length(spots)
            temp(spots(j))=temp(spots(j)-1);
        end
    end
    cp(k)=sum(abs(sign(diff(temp))));
end


if flags.save
    filetit = ['doseopt-paramset-',num2str(flags.paramset),'-ID-',num2str(flags.ID),'-time-',num2str(t_final)];
    save(strcat(pwd,'\Runs\',filetit,'.mat'),'params','dose_times','peak','peak_time','tot_imm','delT','A','Time')
end
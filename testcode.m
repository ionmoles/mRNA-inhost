%% PLOTTING STUFF
% close all
% clear variables
% 
% flags.paramset = 0; %choose the parameters to run the model
% flags.plot = 1; %choose to plot (1) or not (0)
% flags.doses = 1; %number of doses
% 
% run_params;
% 
% scales = [1;params.V0;params.T0;params.B0;params.A0;params.C0;params.F0;params.I0];
% 
% t_final = 256;
% 
% y0 = [1;0;0;0;params.Ai;0;params.Fi;params.Ii];
% 
% odeopts = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-9);
% 
% flags.modtype = 'scaled';
% [ts,ys] = ode23t(@ODEf,[0,t_final/params.t0],y0./scales,odeopts,params,flags);
% 
% params.alphaI = 10*params.alphaI;
% [ts2,ys2] = ode23t(@ODEf,[0,t_final/params.t0],y0./scales,odeopts,params,flags);
% 
% params.alphaI = 100*params.alphaI;
% [ts3,ys3] = ode23t(@ODEf,[0,t_final/params.t0],y0./scales,odeopts,params,flags);
% 
% params.alphaI = 0;
% [ts4,ys4] = ode23t(@ODEf,[0,t_final/params.t0],y0./scales,odeopts,params,flags);
% 
% figure
% plot(ts*params.t0,log10(params.A0*ys(:,5)),'linewidth',2);
% hold on
% plot(ts2*params.t0,log10(params.A0*ys2(:,5)),'linewidth',2);
% plot(ts3*params.t0,log10(params.A0*ys3(:,5)),'linewidth',2);
% plot(ts4*params.t0,log10(params.A0*ys4(:,5)),'linewidth',2);
% DATA = load('Data.mat');
% plot(DATA.Xp,DATA.Yp,'kv');
% plot(DATA.Xm,DATA.Ym,'gs');
% plot(DATA.Xb,DATA.Yb,'r^');
% hold off
% legend('\gamma_I','10\gamma_I','100\gamma_I','\gamma_I=0','location','best');
% 
% figure
% plot(ts*params.t0,params.I0*ys(:,8),'linewidth',2);
% hold on
% plot(ts2*params.t0,params.I0*ys2(:,8),'linewidth',2);
% plot(ts3*params.t0,params.I0*ys3(:,8),'linewidth',2);
% plot(ts4*params.t0,params.I0*ys4(:,8),'linewidth',2);
% hold off
% legend('\gamma_I','10\gamma_I','100\gamma_I','\gamma_I=0','location','best');

%% 
close all
clear variables
flags.paramset = 100;

flags.doses = 1; %number of doses
flags.save = 0; %to save (1) or not (0)
%% 

% alphaB = [12.9/100,12.9/10,12.9];%linspace(1E-2,1E2,100);
lamb = [1.24E-2,80];
I0 = [4.07,0];

dataID = 8;

figure
hold on
for J = 1:length(lamb)
    
    if mod(J,10)==0
        fprintf('Running alpha index %i\n',J);
    end
    
    flags.ID = dataID;
    run_params;
    
    params.lambdaB = lamb(J);
    params.Ii = I0(J);

    scales = [1;params.V0;params.T0;params.B0;params.A0;params.C0;params.F0;params.I0];

    t_final = 256;

    y0 = [1;0;0;0;params.Ai;0;params.Fi;params.Ii];

    odeopts = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-9);

    flags.modtype = 'scaled';

    params.alphaL = 0;
    params.lambdaC=0;
    params.lambdaF=0;
    params.kappaI=0;

    [ts,ys] = ode23t(@ODEf,[0,t_final/params.t0],y0./scales,odeopts,params,flags);
    
    Ts{J} = ts;
    Is{J} = ys(:,8);
    
    plot(ts*params.t0,ys(:,8)*params.I0,'linewidth',2);
    
end

realDATA = load(strcat(pwd,'\data\DATA_ID-',num2str(dataID),'.mat'));
plot(realDATA.X,realDATA.Y,'o');

hold off
xlim([0,30])
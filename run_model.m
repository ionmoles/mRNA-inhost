close all
clear variables

flags.paramset = 100; %choose the parameters to run the model
%paramset
% -1 - artificial parameters for non-dimensional model
% 0 - population parameters from Chapin model
% 1 - Pfizer parameters
% 2 - Moderna parameters
% 100 - these are the different data sets from the literature given by

flags.ID = 8;
flags.plot = 1; %choose to plot (1) or not (0)
flags.datacomp = 1; %choose to plot data comparison (1) or not (0)
flags.doses = 1; %number of doses
flags.save = 0; %to save (1) or not (0)

dose_times = [23;23;28;22;20;22;20;22;22;22;28;28;28;22;28;28;28;28;28;28];
%% 

run_params;

scales = [1;params.V0;params.T0;params.B0;params.A0;params.C0;params.F0;params.I0];

t_final = 256;%1500;%256;
t_dose = dose_times(flags.ID);

y0 = [1;0;0;0;params.Ai;0;params.Fi;params.Ii];

odeopts = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-9);

if flags.paramset==0
    flags.modtype = 'full';
else
    flags.modtype = 'scaled';
end

switch flags.doses
    case 1
        [t,y] = ode23t(@ODEf,[0,t_final],y0,odeopts,params,flags);
    case 2
        sol = ode23t(@ODEf,[0,t_dose],y0,odeopts,params,flags);
        y02 = sol.y(:,end);
        y02(1) = 1;
        sol = odextend(sol,[],t_final,y02,odeopts,params,flags);
        t = sol.x';
        y = sol.y';
end

L = y(:,1);
V = y(:,2);
T = y(:,3);
B = y(:,4);
A = y(:,5);
C = y(:,6);
F = y(:,7);
I = y(:,8);

flags.modtype = 'scaled';

params.alphaL = 0;
params.lambdaC=0;
params.lambdaF=0;
params.kappaI=0;

switch flags.doses
    case 1
        [ts,ys] = ode23t(@ODEf,[0,t_final/params.t0],y0./scales,odeopts,params,flags);
    case 2
        sol = ode23t(@ODEf,[0,t_dose/params.t0],y0./scales,odeopts,params,flags);
        y02 = sol.y(:,end);
        y02(1) = 1;
        sol = odextend(sol,[],t_final/params.t0,y02,odeopts,params,flags);
        ts = sol.x';
        ys = sol.y';
        dspot = find(abs(ts-t_dose/params.t0)==min(abs(ts-t_dose/params.t0)),1);
        Z = zeros(dspot,1);
end

Ls = ys(:,1);
Vs = ys(:,2);
Ts = ys(:,3);
Bs = ys(:,4);
As = ys(:,5);
Cs = ys(:,6);
Fs = ys(:,7);
Is = ys(:,8);

asdat = get_asymptotics(ts,flags,params,t_final);

if flags.plot

    for k=1:8
        figure
    %     plot(t,y(:,k),'linewidth',2);
        plot(ts,ys(:,k),'linewidth',2);
        hold on
    %     plot(ts*params.t0,ys(:,k)*scales(k),'r--','linewidth',2);
    %     plot(ts2,ys2(:,k),'r--','linewidth',2);
    %     plot(ts,ys(:,k),'r--','linewidth',2);

        switch k
            case 2
    %             plot(ts,V0+params.epsilon*V1,'k-.','linewidth',2);
            case 3
                plot(ts,asdat.T0+params.epsilon*asdat.T1,'m-.','linewidth',2);
                plot(ts,params.epsilon^(-1/3)*(asdat.TI0+params.epsilon^(1/3)*asdat.TI1),'r-.','linewidth',2);
                plot(ts,1/params.epsilon*(asdat.T0hat+params.epsilon*asdat.T1hat),'k-.','linewidth',2);
    %             switch flags.doses
    %                 case 1
    %                     plot(ts,1/params.epsilon*(T0hat+params.epsilon*T1hat),'r--','linewidth',2);
    %                 case 2
    %                     plot(ts,1/params.epsilon*(T0hat+params.epsilon*T1hat)+[Z;T0(1:end-dspot)+params.epsilon*T1(1:end-dspot)],'r--','linewidth',2);
    %             end
    %             plot(ts,params.epsilon*Ttil,'g-.','linewidth',2);
                ylim([0,1.05*max(Ts)]);
            case 4
                plot(ts,asdat.B0+params.epsilon*asdat.B1,'m-.','linewidth',2);
                plot(ts,params.epsilon^(-2/3)*(asdat.BI0+params.epsilon^(1/3)*asdat.BI1),'r-.','linewidth',2);
                plot(ts,1/params.epsilon^2*(asdat.B0hat+params.epsilon*asdat.B1hat),'k-.','linewidth',2);
    %             plot(ts,1/params.epsilon*Btil,'g-.','linewidth',2);
                ylim([0,1.05*max(Bs)]);
    %             ylim([0,0.15/params.epsilon^3]);
            case 5
                plot(ts,asdat.A0+params.epsilon*asdat.A1,'m-.','linewidth',2);
                plot(ts,asdat.AI0+params.epsilon^(1/3)*asdat.AI1,'r-.','linewidth',2);
                plot(ts,1/params.epsilon^2*(asdat.A0hat+params.epsilon*asdat.A1hat),'k-.','linewidth',2);
                ylim([0,1.05*max(As)]);
            case 6
                plot(ts,asdat.C0+params.epsilon*asdat.C1,'m-.','linewidth',2);
                plot(ts,params.epsilon^(-1/3)*(asdat.CI0+params.epsilon^(1/3)*asdat.CI1),'r-.','linewidth',2);
                plot(ts,1/params.epsilon*(asdat.C0hat+params.epsilon*asdat.C1hat),'k-.','linewidth',2);
                ylim([0,1.05*max(Cs)]);
            case 7
                plot(ts,Ts,'r--','linewidth',2);
                ylim([0,1.05*max(Ts)]);
            case 8
                plot(ts,asdat.I0+params.epsilon*asdat.I1,'m-.','linewidth',2);
                plot(ts,params.epsilon^(-2/3)*(asdat.II0+params.epsilon^(1/3)*asdat.II1),'r-.','linewidth',2);
                plot(ts,asdat.I0hat+params.epsilon*asdat.I1hat,'k-.','linewidth',2);
                plot(ts,asdat.Itil,'g-.','linewidth',2);
                ylim([0,1.05*max(Is)]);
        end
        hold off
    end

    if flags.datacomp
        switch flags.paramset
            case 100
                figure
                DATA = load(strcat(pwd,'\data\DATA_ID-',num2str(flags.ID),'.mat'));
                switch flags.ID
                    case {1,2,3,4,5,11,12,13,15,16,17,18,19,20}
                        plot(ts*params.t0,log10(params.A0*As),'linewidth',2);
                        hold on
                        plot(ts*params.t0,log10(params.A0*(asdat.A0+params.epsilon*asdat.A1)),'g-.','linewidth',2);
                        plot(ts*params.t0,log10(params.A0*(asdat.AI0+params.epsilon^(1/3)*asdat.AI1)),'r--','linewidth',2);
                        plot(ts*params.t0,log10(params.A0/params.epsilon^2*(asdat.A0hat+params.epsilon*asdat.A1hat)),'k-.','linewidth',2);
                        ymin = -1;
                        ymax = 8;
                    case {6,7}
                        plot(ts*params.t0,params.F0*Fs,'linewidth',2);
                        hold on
                        plot(ts*params.t0,params.F0*Ts,'r--','linewidth',2);
                        ymin = 0;
                        ymax = 1.1*max(DATA.Y);
                    case {8,9,10,14}
                        plot(ts*params.t0,params.I0*Is,'linewidth',2);
                        hold on
                        plot(ts*params.t0,params.I0*(asdat.I0+params.epsilon*asdat.I1),'m-.','linewidth',2);
                        plot(ts*params.t0,params.I0*params.epsilon^(-2/3)*(asdat.II0+params.epsilon^(1/3)*asdat.II1),'r-.','linewidth',2);
                        plot(ts*params.t0,params.I0*(asdat.I0hat+params.epsilon*asdat.I1hat),'k-.','linewidth',2);
                        plot(ts*params.t0,params.I0*asdat.Itil,'g-.','linewidth',2);
                        ymin = 0;
                        ymax = 1.1*max(DATA.Y);
                end
                
                plot(DATA.X,DATA.Y,'ro');
                hold off
                ylim([ymin,ymax]);
                
            case 0
                figure
                plot(t,log10(A),'linewidth',2);
                hold on
                plot(ts*params.t0,log10(params.A0*As),'m--','linewidth',2);
                plot(ts*params.t0,log10(params.A0*(asdat.A0+params.epsilon*asdat.A1)),'g-.','linewidth',2);
                plot(ts*params.t0,log10(params.A0*(asdat.AI0+params.epsilon^(1/3)*asdat.AI1)),'r--','linewidth',2);
                plot(ts*params.t0,log10(params.A0/params.epsilon^2*(asdat.A0hat+params.epsilon*asdat.A1hat)),'k-.','linewidth',2);
                DATA = load('Data.mat');
                plot(DATA.Xp,DATA.Yp,'kv');
                plot(DATA.Xm,DATA.Ym,'gs');
                plot(DATA.Xb,DATA.Yb,'r^');
                hold off
                ylim([-1,8]);

                real_time = ts*params.t0;
                delspot = find(abs(real_time-params.delta*params.t0)==min(abs(real_time-params.delta*params.t0)));
                onespot = find(abs(real_time-params.t0)==min(abs(real_time-params.t0)));
                intspot = find(abs(real_time-params.epsilon^(-1/3)*params.t0)==min(abs(real_time-params.epsilon^(-1/3)*params.t0)));
                longspot = find(abs(real_time-params.epsilon^(-1)*params.t0)==min(abs(real_time-params.epsilon^(-1)*params.t0)));
                decayspot = find(abs(real_time-(-2/params.alphaT)*params.epsilon^(-1)*log(params.epsilon)*params.t0)==min(abs(real_time-(-2/params.alphaT)*params.epsilon^(-1)*log(params.epsilon)*params.t0)));

                figure
                plot(t,log10(A),'linewidth',2);
                hold on
                plot(DATA.Xp(DATA.Xp<params.delta*params.t0),DATA.Yp(DATA.Xp<params.delta*params.t0),'gv');
                plot(DATA.Xm(DATA.Xm<params.delta*params.t0),DATA.Ym(DATA.Xm<params.delta*params.t0),'gs');
                plot(DATA.Xb(DATA.Xb<params.delta*params.t0),DATA.Yb(DATA.Xb<params.delta*params.t0),'g^');
                plot(real_time(1:delspot),log10(params.Ai)*ones(delspot,1),'g--','linewidth',2);

                plot(DATA.Xp(DATA.Xp>=params.delta*params.t0 & DATA.Xp<params.epsilon^(-1/3)*params.t0),DATA.Yp(DATA.Xp>=params.delta*params.t0 & DATA.Xp<params.epsilon^(-1/3)*params.t0),'mv');
                plot(DATA.Xm(DATA.Xm>=params.delta*params.t0 & DATA.Xm<params.epsilon^(-1/3)*params.t0),DATA.Ym(DATA.Xm>=params.delta*params.t0 & DATA.Xm<params.epsilon^(-1/3)*params.t0),'ms');
                plot(DATA.Xb(DATA.Xb>=params.delta*params.t0 & DATA.Xb<params.epsilon^(-1/3)*params.t0),DATA.Yb(DATA.Xb>=params.delta*params.t0 & DATA.Xb<params.epsilon^(-1/3)*params.t0),'m^');
                plot(real_time(delspot:intspot),log10(params.A0*(asdat.A0(delspot:intspot)+params.epsilon*asdat.A1(delspot:intspot))),'m--','linewidth',2);

                plot(DATA.Xp(DATA.Xp>=params.epsilon^(-1/3)*params.t0 & DATA.Xp<params.epsilon^(-1)*params.t0),DATA.Yp(DATA.Xp>=params.epsilon^(-1/3)*params.t0 & DATA.Xp<params.epsilon^(-1)*params.t0),'rv');
                plot(DATA.Xm(DATA.Xm>=params.epsilon^(-1/3)*params.t0 & DATA.Xm<params.epsilon^(-1)*params.t0),DATA.Ym(DATA.Xm>=params.epsilon^(-1/3)*params.t0 & DATA.Xm<params.epsilon^(-1)*params.t0),'rs');
                plot(DATA.Xb(DATA.Xb>=params.epsilon^(-1/3)*params.t0 & DATA.Xb<params.epsilon^(-1)*params.t0),DATA.Yb(DATA.Xb>=params.epsilon^(-1/3)*params.t0 & DATA.Xb<params.epsilon^(-1)*params.t0),'r^');
                plot(real_time(intspot:longspot),log10(params.A0*(asdat.AI0(intspot:longspot)+params.epsilon^(1/3)*asdat.AI1(intspot:longspot))),'r--','linewidth',2);

                plot(DATA.Xp(DATA.Xp>=params.epsilon^(-1)*params.t0 & DATA.Xp<(-2/params.alphaT)*params.epsilon^(-1)*log(params.epsilon)*params.t0),DATA.Yp(DATA.Xp>=params.epsilon^(-1)*params.t0 & DATA.Xp<(-2/params.alphaT)*params.epsilon^(-1)*log(params.epsilon)*params.t0),'kv');
                plot(DATA.Xm(DATA.Xm>=params.epsilon^(-1)*params.t0 & DATA.Xm<(-2/params.alphaT)*params.epsilon^(-1)*log(params.epsilon)*params.t0),DATA.Ym(DATA.Xm>=params.epsilon^(-1)*params.t0 & DATA.Xm<(-2/params.alphaT)*params.epsilon^(-1)*log(params.epsilon)*params.t0),'ks');
                plot(DATA.Xb(DATA.Xb>=params.epsilon^(-1)*params.t0 & DATA.Xb<(-2/params.alphaT)*params.epsilon^(-1)*log(params.epsilon)*params.t0),DATA.Yb(DATA.Xb>=params.epsilon^(-1)*params.t0 & DATA.Xb<(-2/params.alphaT)*params.epsilon^(-1)*log(params.epsilon)*params.t0),'k^');
                plot(real_time(longspot:decayspot),log10(params.A0/params.epsilon^2*(asdat.A0hat(longspot:decayspot)+params.epsilon*asdat.A1hat(longspot:decayspot))),'k--','linewidth',2);

                plot(DATA.Xp(DATA.Xp>=(-2/params.alphaT)*params.epsilon^(-1)*log(params.epsilon)*params.t0),DATA.Yp(DATA.Xp>=(-2/params.alphaT)*params.epsilon^(-1)*log(params.epsilon)*params.t0),'cv');
                plot(DATA.Xm(DATA.Xm>=(-2/params.alphaT)*params.epsilon^(-1)*log(params.epsilon)*params.t0),DATA.Ym(DATA.Xm>=(-2/params.alphaT)*params.epsilon^(-1)*log(params.epsilon)*params.t0),'cs');
                plot(DATA.Xb(DATA.Xb>=(-2/params.alphaT)*params.epsilon^(-1)*log(params.epsilon)*params.t0),DATA.Yb(DATA.Xb>=(-2/params.alphaT)*params.epsilon^(-1)*log(params.epsilon)*params.t0),'c^');
                plot(real_time(decayspot:end),log10(params.A0/params.epsilon^2*(asdat.A0hat(decayspot:end)+params.epsilon*asdat.A1hat(decayspot:end))),'c--','linewidth',2);
                hold off
        end
    end
end

if flags.save
    if flags.paramset<0
        filetit = ['modeldat-paramset-m',num2str(abs(flags.paramset)),'-ID-',num2str(flags.ID),'-time-',num2str(t_final),'-doses-',num2str(flags.doses)];
    else
        filetit = ['modeldat-paramset-',num2str(flags.paramset),'-ID-',num2str(flags.ID),'-time-',num2str(t_final),'-doses-',num2str(flags.doses)];
    end
    save(strcat(pwd,'\Runs\',filetit,'.mat'),'t','y','ts','ys','params','flags','t_final')
end
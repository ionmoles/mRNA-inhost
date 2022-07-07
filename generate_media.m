close all
clear variables
%% BASIC SETUP

pcols = {1/255*[0,196,255],1/255*[225,12,62],1/255*[55,176,116],1/255*[244,201,107],1/255*[149,52,235],1/255*[240,195,239],1/255*[176,176,176]};
bcol = 1/255*[224,224,224];

curpath = pwd;
plotpath = strcat(curpath,'\plots\');

pos = [2,4.5,6,4.5];

%% 
%Experiment
%-1 - Chapin appendix figures
%0 - comparison of full model to reduced model for artificial parameters
%1 - asymptotic regions compared to full model
%2 - 1 dose model comparison with asymptotics and data
%3 - 1 and 2 dose model comparison with data
%4 - optimization graphs
experiment = 3;

switch(experiment)
    case -1
        DATA = load('Chapin_dat.mat');
        ya = DATA.ya;
        y = DATA.y;
        t = DATA.t;
        
        pnames={'Lapp';'Vapp';'Tapp';'Bapp';'Aapp';'Capp';'Fapp';'Iapp'};
        ynames={'$L$';'$V$';'$T$';'$B$';'$A$';'$C$';'$F$';'$I$'};
        
        fig_num = 0;
        for k=[1,2,3,4,5,6,7,8]
            fig_num = fig_num + 1;
            [fig(fig_num),figprop(fig_num)]=new_figure(pnames{fig_num},pos);
            plot(t,y(:,k),'linewidth',2);
            hold on
            plot(t,ya(:,k),'linestyle','--','linewidth',2);
            hold off
            figprop(fig_num).xlab = xlabel('$t$ (days)');
            figprop(fig_num).ylab = ylabel(ynames{fig_num});
            figprop(fig_num).leg = legend('Original Model','Reduced Model','location','best');
        end
        
        
    case 0
        flags.paramset = -1;

        run_params;

        scales = [1;params.V0;params.T0;params.B0;params.A0;params.C0;params.F0;params.I0];

        t_final = 1500;%256;

        y0 = [1;0;0;0;params.Ai;0;params.Fi;params.Ii];

        odeopts = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-9);
        
        flags.modtype = 'scaled';
        [t,y] = ode23t(@ODEf,[0,t_final],y0,odeopts,params,flags);
        
        flags.modtype = 'scaled';

        params.alphaL = 0;
        params.lambdaC =0;
        params.lambdaF = 0;
       
        [ts,ys] = ode23t(@ODEf,t,y0./scales,odeopts,params,flags);
        
        
        params.alphaL = 0;
        params.lambdaC =0;
        params.lambdaF = 0;
        params.kappaI = 0;
        params.kappaF = 0;

        [ts2,ys2] = ode23t(@ODEf,t,y0./scales,odeopts,params,flags);
        
        pnames = {'Lcomp';'Vcomp';'Tcomp';'Bcomp';'Acomp';'Ccomp';'Fcomp';'Icomp'};
        
        for k=1:length(pnames)
            [fig(k),figprop(k)]=new_figure(pnames{k},pos);
            plot(t,y(:,k),'Color',pcols{1},'linewidth',2);
            hold on
            plot(ts*params.t0,ys(:,k)*scales(k),'LineStyle','--','Color',pcols{2},'linewidth',2);
            plot(ts2*params.t0,ys2(:,k)*scales(k),'LineStyle',':','Color',pcols{3},'linewidth',2);
            hold off
            figprop(k).leg = legend({'Full','Reduced','Reduced $\kappa_{\rm I}=0$'},'location','best');
%             figprop(k).xlab = xlabel('$t$ (days)');
            figprop(k).xlab = xlabel('$t$');
            figprop(k).ax = gca;
            if k==1
                xlim([0,50]);
            elseif k==8
                xlim([0,1500]);
            else
                xlim([0,1000]);
            end
        end
    case 1
        DATA = load(strcat(pwd,'\Runs\','modeldat-paramset-m1-time-1500-doses-1'));
        params = DATA.params;
        flags = DATA.flags;
        t_final = DATA.t_final;
        t = DATA.t;
        y = DATA.y;
        ys = DATA.ys;
        
        L = DATA.y(:,1);
        V = DATA.y(:,2);
        T = DATA.y(:,3);
        B = DATA.y(:,4);
        A = DATA.y(:,5);
        C = DATA.y(:,6);
        F = DATA.y(:,7);
        I = DATA.y(:,8);
        
        ts = DATA.ts;
        Ls = DATA.ys(:,1);
        Vs = DATA.ys(:,2);
        Ts = DATA.ys(:,3);
        Bs = DATA.ys(:,4);
        As = DATA.ys(:,5);
        Cs = DATA.ys(:,6);
        Fs = DATA.ys(:,7);
        Is = DATA.ys(:,8);
        
        t1 = params.delta;
        t2 = 1;
        t3 = params.epsilon^(-1/3);
        t4 = 1/params.epsilon;
        t5 = -2/params.epsilon/params.alphai*log(params.epsilon);
        
        pnames={'Lass';'Vass';'Tass';'Bass';'Aass';'Cass';'Fass';'Iass'};
        
        sf = 1.1;
        
        asdat=get_asymptotics(ts,flags,params,t_final);
        ya2(:,1) = exp(-ts);
        ya2(:,2) = (exp(-params.epsilon*params.alphaV*ts)-exp(-ts))/(1-params.epsilon*params.alphaV);
        ya2(:,3) = exp(-params.epsilon*params.alphaT*ts)/(1-params.epsilon*params.alphaV).*((exp(-(1-params.epsilon*params.alphaT)*ts)-1)/(1-params.epsilon*params.alphaT)-(exp(-params.epsilon*(params.alphaV-params.alphaT)*ts)-1)/(params.epsilon*(params.alphaV-params.alphaT)));
        ya2(:,4) = asdat.B0+params.epsilon*asdat.B1;
        ya2(:,5) = asdat.A0+params.epsilon*asdat.A1;
        ya2(:,6) = exp(-params.epsilon*params.alphaC*ts)/(1-params.epsilon*params.alphaV).*((exp(-(1-params.epsilon*params.alphaC)*ts)-1)/(1-params.epsilon*params.alphaC)-(exp(-params.epsilon*(params.alphaV-params.alphaC)*ts)-1)/(params.epsilon*(params.alphaV-params.alphaC)));
        ya2(:,7) = ya2(:,3);
        ya2(:,8) = asdat.I0+params.epsilon*asdat.I1;
        
        ya3(:,1) = zeros(length(ts),1);
        ya3(:,2) = zeros(length(ts),1);
        ya3(:,3) = zeros(length(ts),1);
        ya3(:,4) = params.epsilon^(-2/3)*(asdat.BI0+0*params.epsilon^(1/3)*asdat.BI1);
        ya3(:,5) = asdat.AI0+0*params.epsilon^(1/3)*asdat.AI1;
        ya3(:,6) = zeros(length(ts),1);
        ya3(:,7) = zeros(length(ts),1);
        ya3(:,8) = params.epsilon^(-2/3)*(asdat.II0+0*params.epsilon^(1/3)*asdat.II1);
        
        ya4(:,1) = zeros(length(ts),1);
        ya4(:,2) = zeros(length(ts),1);
        ya4(:,3) = zeros(length(ts),1);
        ya4(:,4) = 1/params.epsilon^2*(asdat.B0hat+0*params.epsilon*asdat.B1hat);
        ya4(:,5) = 1/params.epsilon^2*(asdat.A0hat+0*params.epsilon*asdat.A1hat);
        ya4(:,6) = zeros(length(ts),1);
        ya4(:,7) = zeros(length(ts),1);
        ya4(:,8) = asdat.I0hat+0*params.epsilon*asdat.I1hat;
        
        ya5(:,1) = zeros(length(ts),1);
        ya5(:,2) = zeros(length(ts),1);
        ya5(:,3) = zeros(length(ts),1);
        ya5(:,4) = zeros(length(ts),1);
        ya5(:,5) = zeros(length(ts),1);
        ya5(:,6) = zeros(length(ts),1);
        ya5(:,7) = zeros(length(ts),1);
        ya5(:,8) = asdat.Itil;
        
        plabs = {'$L$';'$V$';'$T$';'$B$';'$A$';'$C$';'$F$';'$I$'};
        legplots = {[1,5,6];[1,5,6];[1,5,6];[1,5,6,7,8];[1,5,6,7,8];[1,5,6];[1,5,6];[1,5,6,7,8,9]};
        legentries = {{'Numeric (full)','Numeric (reduced)','Analytic'},{'Numeric (full)','Numeric (reduced)','Analytic'},{'Numeric (full)','Numeric (reduced)','Analytic'},{'Numeric (full)','Numeric (reduced)','$t_2$','$t_3$','$t_4$'},{'Numeric (full)','Numeric (reduced)','$t_2$','$t_3$','$t_4$'},{'Numeric (full)','Numeric (reduced)','Analytic'},{'Numeric (full)','Numeric (reduced)','Analytic'},{'Numeric (full)','Numeric (reduced)','$t_2$','$t_3$','$t_4$','$t_5$'}};
        
        x_limit = [t3,1.5*t4,1.5*t4,1.5*t4,2.5*t4,1.5*t4,1.5*t4,0.5*t4];
        
        fignum = 0;
        for k=1:length(pnames)
            p = 0;
            fignum = fignum + 1;
            [fig(k),figprop(k)]=new_figure(pnames{k},pos);
            p(1)=plot(t,y(:,k),'Color',pcols{6},'linewidth',2);
            UB = max(y(:,k))*sf;
            hold on
            p(2)=patch([t1,t2,t2,t1],[UB,UB,0,0],bcol);
            p(3)=patch([t3,t4,t4,t3],[UB,UB,0,0],bcol);
            p(4)=patch([t5,ts(end),ts(end),t5],[UB,UB,0,0],bcol);
            
            p(1)=plot(t,y(:,k),'Color',pcols{6},'linewidth',2);

            p(5)=plot(ts,ys(:,k),'Color',pcols{1},'linewidth',2);
            
            if ismember(6,legplots{k})
                p(6)=plot(ts,ya2(:,k),'LineStyle','--','Color',pcols{2},'linewidth',2);
            end
            if ismember(7,legplots{k})
                p(7)=plot(ts,ya3(:,k),'LineStyle','--','Color',pcols{3},'linewidth',2);
            end
            if ismember(8,legplots{k})
                p(8)=plot(ts,ya4(:,k),'LineStyle','--','Color',pcols{4},'linewidth',2);
            end
            if ismember(9,legplots{k})
                p(9)=plot(ts,ya5(:,k),'LineStyle','--','Color',pcols{5},'linewidth',2);
            end

            hold off
            xlim([0,x_limit(k)])
            ylim([0,UB])
            figprop(k).xlab=xlabel('$t$');
            figprop(k).ylab=ylabel(plabs{k});
            figprop(k).leg = legend(p(legplots{k}),legentries{k},'location','best');
        end
        
        fignum = fignum + 1;
        [fig(fignum),figprop(fignum)]=new_figure('Iassout',pos);
        p = 0;
        p(1)=plot(t,y(:,8),'Color',pcols{6},'linewidth',2);
        UB = max(y(:,8))*sf;
        hold on
        p(2)=patch([t1,t2,t2,t1],[UB,UB,0,0],bcol);
        p(3)=patch([t3,t4,t4,t3],[UB,UB,0,0],bcol);
        p(4)=patch([t5,ts(end),ts(end),t5],[UB,UB,0,0],bcol);
        
        p(1)=plot(t,y(:,8),'Color',pcols{6},'linewidth',2);

        p(5)=plot(ts,ys(:,8),'Color',pcols{1},'linewidth',2);

        p(6)=plot(ts,ya2(:,8),'LineStyle','--','Color',pcols{2},'linewidth',2);
        p(7)=plot(ts,ya3(:,8),'LineStyle','--','Color',pcols{3},'linewidth',2);
        p(8)=plot(ts,ya4(:,8),'LineStyle','--','Color',pcols{4},'linewidth',2);
        p(9)=plot(ts,ya5(:,8),'LineStyle','--','Color',pcols{5},'linewidth',2);

        hold off
        ylim([0,UB])
        figprop(fignum).xlab=xlabel('$t$');
        figprop(fignum).ylab=ylabel(plabs{8});
        figprop(fignum).leg = legend(p(legplots{8}),legentries{8},'location','best');
        
        fignum = fignum + 1;
        [fig(fignum),figprop(fignum)]=new_figure('Aassin',pos);
        p = 0;
        p(1)=plot(t,y(:,5),'Color',pcols{6},'linewidth',2);
        UB = max(y(:,5))*sf;
        hold on
        p(2)=patch([t1,t2,t2,t1],[UB,UB,0,0],bcol);
        p(3)=patch([t3,t4,t4,t3],[UB,UB,0,0],bcol);
        p(4)=patch([t5,ts(end),ts(end),t5],[UB,UB,0,0],bcol);
        
        p(1)=plot(t,y(:,5),'Color',pcols{6},'linewidth',2);

        p(5)=plot(ts,ys(:,5),'Color',pcols{1},'linewidth',2);

        p(6)=plot(ts,ya2(:,5),'LineStyle','--','Color',pcols{2},'linewidth',2);
        p(7)=plot(ts,ya3(:,5),'LineStyle','--','Color',pcols{3},'linewidth',2);
        p(8)=plot(ts,ya4(:,5),'LineStyle','--','Color',pcols{4},'linewidth',2);

        hold off
        ylim([0,10])
        xlim([0,15]);
        figprop(fignum).xlab=xlabel('$t$');
        figprop(fignum).ylab=ylabel(plabs{5});
        figprop(fignum).leg = legend(p(legplots{5}),legentries{5},'location','best');
    case 2
        fignum = 0;
        maxDATA = load('maxtimes.mat');
        for k=1:20
            clear ya2 ya3 ya4 ya5
            pnames={'A1';'A2';'A3';'A4';'A5';'F6';'F7';'I8';'I9';'I10';'A11';'A12';'A13';'I14';'A15';'A16';'A17';'A18';'A19';'A20'};
            plotspot = [5;5;5;5;5;7;7;8;8;8;5;5;5;8;5;5;5;5;5;5];
            
            DATA2D = load(strcat(pwd,'\Runs\','modeldat-patientID-',num2str(k),'.mat'));
            DATA = load(strcat(pwd,'\Runs\','modeldat-patientID-',num2str(k),'-1dose.mat'));
            realDATA = load(strcat(pwd,'\data\DATA_ID-',num2str(k),'.mat'));
            
            params = DATA.params;
            flags = DATA.flags;
            t_final = DATA.t_final;
            t_dose = DATA.t_dose;

            ts = DATA.ts;
            ys = DATA.ys;
            
            scales = [1;params.V0;params.T0;params.B0;params.A0;params.C0;params.F0;params.I0];

            t = params.t0*ts;
            
            t1 = params.delta*params.t0;
            t2 = 1*params.t0;
            t3 = params.epsilon^(-1/3)*params.t0;
            t4 = 1/params.epsilon*params.t0;
            t5 = -2/params.epsilon/params.alphai*log(params.epsilon)*params.t0;
            
            t2spot = find(abs(t-t2)==min(abs(t-t2)));
            t3spot = find(abs(t-t3)==min(abs(t-t3)));
            t4spot = find(abs(t-t4)==min(abs(t-t4)));
            if min([params.alphaV,params.alphaT,params.alphaB])~=params.alphaB
                t5spot = find(abs(t-t5)==min(abs(t-t5)));
            else
                t5spot = length(ts);
            end


            sf = 1.1;

            asdat=get_asymptotics(ts,flags,params,t_final);
            ya2(:,1) = exp(-ts);
            ya2(:,2) = (exp(-params.epsilon*params.alphaV*ts)-exp(-ts))/(1-params.epsilon*params.alphaV);
            ya2(:,3) = exp(-params.epsilon*params.alphaT*ts)/(1-params.epsilon*params.alphaV).*((exp(-(1-params.epsilon*params.alphaT)*ts)-1)/(1-params.epsilon*params.alphaT)-(exp(-params.epsilon*(params.alphaV-params.alphaT)*ts)-1)/(params.epsilon*(params.alphaV-params.alphaT)));
            ya2(:,4) = asdat.B0+params.epsilon*asdat.B1;
            ya2(:,5) = asdat.A0+params.epsilon*asdat.A1;
            ya2(:,6) = exp(-params.epsilon*params.alphaC*ts)/(1-params.epsilon*params.alphaV).*((exp(-(1-params.epsilon*params.alphaC)*ts)-1)/(1-params.epsilon*params.alphaC)-(exp(-params.epsilon*(params.alphaV-params.alphaC)*ts)-1)/(params.epsilon*(params.alphaV-params.alphaC)));
            ya2(:,7) = ya2(:,3);
            ya2(:,8) = asdat.I0+params.epsilon*asdat.I1;

            ya3(:,1) = zeros(length(ts),1);
            ya3(:,2) = zeros(length(ts),1);
            ya3(:,3) = zeros(length(ts),1);
            ya3(:,4) = params.epsilon^(-2/3)*(asdat.BI0+0*params.epsilon^(1/3)*asdat.BI1);
            ya3(:,5) = asdat.AI0+0*params.epsilon^(1/3)*asdat.AI1;
            ya3(:,6) = zeros(length(ts),1);
            ya3(:,7) = zeros(length(ts),1);
            ya3(:,8) = params.epsilon^(-2/3)*(asdat.II0+0*params.epsilon^(1/3)*asdat.II1);

            ya4(:,1) = zeros(length(ts),1);
            ya4(:,2) = zeros(length(ts),1);
            ya4(:,3) = zeros(length(ts),1);
            ya4(:,4) = 1/params.epsilon^2*(asdat.B0hat+0*params.epsilon*asdat.B1hat);
            ya4(:,5) = 1/params.epsilon^2*(asdat.A0hat+0*params.epsilon*asdat.A1hat);
            ya4(:,6) = zeros(length(ts),1);
            ya4(:,7) = zeros(length(ts),1);
            ya4(:,8) = asdat.I0hat+0*params.epsilon*asdat.I1hat;

            ya5(:,1) = zeros(length(ts),1);
            ya5(:,2) = zeros(length(ts),1);
            ya5(:,3) = zeros(length(ts),1);
            ya5(:,4) = zeros(length(ts),1);
            ya5(:,5) = zeros(length(ts),1);
            ya5(:,6) = zeros(length(ts),1);
            ya5(:,7) = zeros(length(ts),1);
            ya5(:,8) = asdat.Itil;

            plabs={'$\log{A}$';'$\log{A}$';'$\log{A}$';'$\log{A}$';'$\log{A}$';'$\log{F}$';'$\log{F}$';'$I$';'$I$';'$I$';'$\log{A}$';'$\log{A}$';'$\log{A}$';'$I$';'$\log{A}$';'$\log{A}$';'$\log{A}$';'$\log{A}$';'$\log{A}$';'$\log{A}$'};
            if min([params.alphaV,params.alphaT,params.alphaB])~=params.alphaB
                legplots{k} =[1,5,6,7,8,9];
                legentries{k} = {'Model','Data','$t_2$','$t_3$','$t_4$','$t_5$'};
            else
                legplots{k}=[1,5,6,7,8];
                legentries{k} = {'Model','Data','$t_2$','$t_3$','$t_4$'};
            end

            x_limit = max(realDATA.X)*ones(1,20);
%             x_limit(6) = 22;
%             x_limit(8) = 22;
%             x_limit(9) = 22;
%             x_limit(10) = 22;
            
            p = 0;
            fignum = fignum + 1;
            [fig(fignum),figprop(fignum)]=new_figure(pnames{k},pos);
            if plotspot(k)==5 || plotspot(k)==7
                p(1)=plot(params.t0*ts,log10(scales(plotspot(k))*ys(:,plotspot(k))),'Color',pcols{1},'linewidth',2);
                UB = max(max(log10(scales(plotspot(k))*ys(:,plotspot(k)))),max(realDATA.Y))*sf;
                LB = -1;
            else
                p(1)=plot(params.t0*ts,scales(plotspot(k))*ys(:,plotspot(k)),'Color',pcols{1},'linewidth',2);
                UB = max(max(scales(plotspot(k))*ys(:,plotspot(k))),max(realDATA.Y))*sf;
                LB = 0;
            end
            
            hold on
            p(2)=patch([t1,t2,t2,t1],[UB,UB,LB,LB],bcol);
            p(3)=patch([t3,t4,t4,t3],[UB,UB,LB,LB],bcol);
            if min([params.alphaV,params.alphaT,params.alphaB])~=params.alphaB
                p(4)=patch([t5,params.t0*ts(end),params.t0*ts(end),t5],[UB,UB,LB,LB],bcol);
            end

            if plotspot(k)==5 
                p(1)=plot(params.t0*ts,log10(scales(plotspot(k))*ys(:,plotspot(k))),'Color',pcols{1},'linewidth',2);
                p(5) = plot(realDATA.X,realDATA.Y,'ko','MarkerSize',10);
                p(6)=plot(params.t0*ts(t2spot:t3spot),log10(scales(plotspot(k))*ya2(t2spot:t3spot,plotspot(k))),'LineStyle','--','Color',pcols{2},'linewidth',2);
                p(7)=plot(params.t0*ts(t3spot:t4spot),log10(scales(plotspot(k))*ya3(t3spot:t4spot,plotspot(k))),'LineStyle','--','Color',pcols{3},'linewidth',2);
                p(8)=plot(params.t0*ts(t3spot:t5spot),log10(scales(plotspot(k))*ya4(t3spot:t5spot,plotspot(k))),'LineStyle','--','Color',pcols{4},'linewidth',2);
            elseif plotspot(k)==7
                p(1)=plot(params.t0*ts,log10(scales(plotspot(k))*ys(:,plotspot(k))),'Color',pcols{1},'linewidth',2);
                p(5) = plot(realDATA.X,realDATA.Y,'ko','MarkerSize',10);
                p(6)=plot(params.t0*ts(t2spot:t3spot),log10(scales(plotspot(k))*ya2(t2spot:t3spot,plotspot(k))),'LineStyle','--','Color',pcols{2},'linewidth',2);
                p(7)=plot(params.t0*ts(t3spot:t4spot),log10(scales(plotspot(k))*ya2(t3spot:t4spot,plotspot(k))),'LineStyle','--','Color',pcols{3},'linewidth',2);
                p(8)=plot(params.t0*ts(t3spot:t5spot),log10(scales(plotspot(k))*ya2(t3spot:t5spot,plotspot(k))),'LineStyle','--','Color',pcols{4},'linewidth',2);
            else
                p(1)=plot(params.t0*ts,scales(plotspot(k))*ys(:,plotspot(k)),'Color',pcols{1},'linewidth',2);
                p(5) = plot(realDATA.X,realDATA.Y,'ko','MarkerSize',10);
                p(6)=plot(params.t0*ts(t2spot:t3spot),scales(plotspot(k))*ya2(t2spot:t3spot,plotspot(k)),'LineStyle','--','Color',pcols{2},'linewidth',2);
                p(7)=plot(params.t0*ts(t3spot:t4spot),scales(plotspot(k))*ya3(t3spot:t4spot,plotspot(k)),'LineStyle','--','Color',pcols{3},'linewidth',2);
                p(8)=plot(params.t0*ts(t3spot:t5spot),scales(plotspot(k))*ya4(t3spot:t5spot,plotspot(k)),'LineStyle','--','Color',pcols{4},'linewidth',2);
            end
            

            if min([params.alphaV,params.alphaT,params.alphaB])~=params.alphaB
                if plotspot(k)==8
                    p(9)=plot(params.t0*ts(t5spot:end),scales(plotspot(k))*ya5(t5spot:end,plotspot(k)),'LineStyle','--','Color',pcols{5},'linewidth',2);
                else
                    p(9)=plot(params.t0*ts(t5spot:end),log10(scales(plotspot(k))*ya4(t5spot:end,plotspot(k))),'LineStyle','--','Color',pcols{5},'linewidth',2);
                end
            end
            
            p(10) = plot(maxDATA.tImax(k)*ones(1000,1),linspace(LB,UB,1000),'Linestyle','--','Color',pcols{end-1},'linewidth',2);
            p(11) = plot(maxDATA.tBmax(k)*ones(1000,1),linspace(LB,UB,1000),'Linestyle','--','Color',pcols{end-1},'linewidth',2);
            p(12) = plot(maxDATA.tAmax(k)*ones(1000,1),linspace(LB,UB,1000),'Linestyle','--','Color',pcols{end-1},'linewidth',2);

            hold off
            xlim([0,x_limit(k)])
            ylim([LB,UB])
            figprop(fignum).xlab=xlabel('$t$ (days)');
            figprop(fignum).ylab=ylabel(plabs{k});
            figprop(fignum).leg = legend(p(legplots{k}),legentries{k},'location','best');
        end
        
    case 3
        fignum = 0;
        for k=1:20
            pnames={'2DA1';'2DA2';'2DA3';'2DA4';'2DA5';'2DF6';'2DF7';'2DI8';'2DI9';'2DI10';'2DA11';'2DA12';'2DA13';'2DI14';'2DA15';'2DA16';'2DA17';'2DA18';'2DA19';'2DA20'};
            plotspot = [5;5;5;5;5;7;7;8;8;8;5;5;5;8;5;5;5;5;5;5];
            
            DATA2D = load(strcat(pwd,'\Runs\','modeldat-patientID-',num2str(k),'.mat'));
            DATA = load(strcat(pwd,'\Runs\','modeldat-patientID-',num2str(k),'-1dose.mat'));
            realDATA = load(strcat(pwd,'\data\DATA_ID-',num2str(k),'.mat'));
            params = DATA.params;
            params2D = DATA2D.params;
            
            flags = DATA.flags;
            flags2D = DATA2D.flags;
            
            t_final = DATA.t_final;
            t_dose = DATA.t_dose;
            t_final2D = DATA2D.t_final;
            t_dose2D = DATA2D.t_dose;

            ts = DATA.ts;
            ts2D = DATA2D.ts;
            ys = DATA.ys;
            ys2D = DATA2D.ys;
            
            scales = [1;params.V0;params.T0;params.B0;params.A0;params.C0;params.F0;params.I0];
            scales2D = [1;params2D.V0;params2D.T0;params2D.B0;params2D.A0;params2D.C0;params2D.F0;params2D.I0];

            t = params.t0*ts;
            t2D = params2D.t0*ts2D;
            
            t1 = params.delta*params.t0;
            t2 = 1*params.t0;
            t3 = params.epsilon^(-1/3)*params.t0;
            t4 = 1/params.epsilon*params.t0;
            t5 = -2/params.epsilon/params.alphai*log(params.epsilon)*params.t0;
            
            sf = 1.1;

            
            plabs={'$\log{A}$';'$\log{A}$';'$\log{A}$';'$\log{A}$';'$\log{A}$';'$\log{F}$';'$\log{F}$';'$I$';'$I$';'$I$';'$\log{A}$';'$\log{A}$';'$\log{A}$';'$I$';'$\log{A}$';'$\log{A}$';'$\log{A}$';'$\log{A}$';'$\log{A}$';'$\log{A}$'};
            if k==3
                legentries{k} = {'Model-1 dose','Model-2 dose','Data'};
                legplots{k} =[1,5,6];
            else
                legentries{k} = {'Model-1 dose','Model-2 dose','Data','2nd dose'};
                legplots{k} =[1,5,6,7];
            end

            x_limit = max(realDATA.X)*ones(1,20);
%             x_limit(6) = 22;
%             x_limit(8) = 22;
%             x_limit(9) = 22;
%             x_limit(10) = 22;
            
            p = 0;
            fignum = fignum + 1;
            [fig(fignum),figprop(fignum)]=new_figure(pnames{k},pos);
            if plotspot(k)==5 || plotspot(k)==7
                p(1)=plot(params.t0*ts,log10(scales(plotspot(k))*ys(:,plotspot(k))),'Color',pcols{1},'linewidth',2);
                UB = max(max(log10(scales(plotspot(k))*ys(:,plotspot(k)))),max(realDATA.Y))*sf;
                LB = -1;
            else
                p(1)=plot(params.t0*ts,scales(plotspot(k))*ys(:,plotspot(k)),'Color',pcols{1},'linewidth',2);
                UB = max(max(scales(plotspot(k))*ys(:,plotspot(k))),max(realDATA.Y))*sf;
                LB = 0;
            end
            
            hold on
            p(2)=patch([t1,t2,t2,t1],[UB,UB,LB,LB],bcol);
            p(3)=patch([t3,t4,t4,t3],[UB,UB,LB,LB],bcol);
            if min([params.alphaV,params.alphaT,params.alphaB])~=params.alphaB
                p(4)=patch([t5,params.t0*ts(end),params.t0*ts(end),t5],[UB,UB,LB,LB],bcol);
            end

            if plotspot(k)==5 
                p(1)=plot(params.t0*ts,log10(scales(plotspot(k))*ys(:,plotspot(k))),'Color',pcols{1},'linewidth',2);
                p(5)=plot(params2D.t0*ts2D,log10(scales2D(plotspot(k))*ys2D(:,plotspot(k))),'--','Color',pcols{2},'linewidth',2);
                p(6) = plot(realDATA.X,realDATA.Y,'ko','MarkerSize',10);
            elseif plotspot(k)==7
                p(1)=plot(params.t0*ts,log10(scales(plotspot(k))*ys(:,plotspot(k))),'Color',pcols{1},'linewidth',2);
                p(5)=plot(params2D.t0*ts2D,log10(scales2D(plotspot(k))*ys2D(:,plotspot(k))),'--','Color',pcols{2},'linewidth',2);
                p(6) = plot(realDATA.X,realDATA.Y,'ko','MarkerSize',10);
            else
                p(1)=plot(params.t0*ts,scales(plotspot(k))*ys(:,plotspot(k)),'Color',pcols{1},'linewidth',2);
                p(5)=plot(params2D.t0*ts2D,scales2D(plotspot(k))*ys2D(:,plotspot(k)),'--','Color',pcols{2},'linewidth',2);
                p(6) = plot(realDATA.X,realDATA.Y,'ko','MarkerSize',10);
            end
            p(7) = plot(t_dose2D*ones(1000,1),linspace(LB,UB,1000),'-.','Color',pcols{3},'linewidth',2);
            

            hold off
            xlim([0,x_limit(k)])
            ylim([LB,UB])
            figprop(fignum).xlab=xlabel('$t$ (days)');
            figprop(fignum).ylab=ylabel(plabs{k});
            figprop(fignum).leg = legend(p(legplots{k}),legentries{k},'location','best');
        end
    case 4
        DATA = load(strcat(pwd,'\Runs\','doseopt-paramset-100-ID-3-time-365.mat'));
        params = DATA.params;
        [fig(1),figprop(1)]=new_figure('opt1',pos);
        plot(DATA.dose_times,DATA.delT(2:end)*params.t0,'Color',pcols{1},'linewidth',2);
        xlim([0,350])
        figprop(1).xlab=xlabel('$\tau$ (days)');
        figprop(1).ylab=ylabel('$\Delta t_{\rm FWHM}$');
        
        [fig(2),figprop(2)]=new_figure('opt2',pos);
        plot(DATA.Time{5}*params.t0,DATA.A{5}*params.A0,'Color',pcols{2},'linewidth',2);
        figprop(2).xlab=xlabel('$t$ (days)');
        figprop(2).ylab=ylabel('$A$');
        
        [fig(3),figprop(3)]=new_figure('opt3',pos);
        plot(DATA.Time{50}*params.t0,DATA.A{50}*params.A0,'Color',pcols{2},'linewidth',2);
        figprop(3).xlab=xlabel('$t$ (days)');
        figprop(3).ylab=ylabel('$A$');
        
        [fig(4),figprop(4)]=new_figure('opt4',pos);
        plot(DATA.Time{100}*params.t0,DATA.A{100}*params.A0,'Color',pcols{2},'linewidth',2);
        figprop(4).xlab=xlabel('$t$ (days)');
        figprop(4).ylab=ylabel('$A$');
        
        [fig(5),figprop(5)]=new_figure('opt5',pos);
        plot(DATA.Time{200}*params.t0,DATA.A{200}*params.A0,'Color',pcols{2},'linewidth',2);
        figprop(5).xlab=xlabel('$t$ (days)');
        figprop(5).ylab=ylabel('$A$');
        
end
%% 

 figure_revise(pos,plotpath,fig,figprop);
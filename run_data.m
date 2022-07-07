close all
clear variables
flags.paramset = 100;

flags.doses = 1; %number of doses
flags.save = 1; %to save (1) or not (0)
%% 

dose_times = [23;23;28;22;20;22;20;22;22;22;28;28;28;22;28;28;28;28;28;28];

for datID = 1:20
    
    fprintf('Running model for patient %i\n',datID);
    
    flags.ID = datID;
    run_params;

    scales = [1;params.V0;params.T0;params.B0;params.A0;params.C0;params.F0;params.I0];

    t_final = 256;
    t_dose = dose_times(datID);

    y0 = [1;0;0;0;params.Ai;0;params.Fi;params.Ii];

    odeopts = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-9);

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
    end

    Ls = ys(:,1);
    Vs = ys(:,2);
    Ts = ys(:,3);
    Bs = ys(:,4);
    As = ys(:,5);
    Cs = ys(:,6);
    Fs = ys(:,7);
    Is = ys(:,8);
    
    asdat=get_asymptotics(ts,flags,params,t_final);


    if flags.save
        if flags.doses==1
            filetit = ['modeldat-patientID-',num2str(flags.ID),'-1dose'];
        else
            filetit = ['modeldat-patientID-',num2str(flags.ID)];
        end
        save(strcat(pwd,'\Runs\',filetit,'.mat'),'ts','ys','params','flags','t_final','t_dose')
    end
end
close all
clear variables

experiment = 2;
%0 - produce an excel file for non-dimensional parameters
%1 - produce data counts
%2 - produce excel file for optima

switch experiment
    case 0
        row_tits = {'ID','V0','T0','B0','A0','C0','F0','I0','t0','epsilon','delta','alphaL','alphaV','alphaT','alphaB','alphaC','alphaI','lambdaB','lambdaC','lambdaF','lambdaI','kappaF','kappaI','A','F','I'};
        P = zeros(20,26);
        
        t_tits = {'ID','t1','t2','t3','t4','t5','tau'};
        T_scales = zeros(20,5);
        tau=[23;23;28;22;20;22;20;22;22;22;28;28;28;22;28;28;28;28;28;28];
        %Estimated 23 for first two based on spread of second dose
        %Wang - determined based on measurements between 3 and 14 weeks
        %post second dose
        flags.paramset = 100;
        for k=1:20
            flags.ID = k;
            run_params;
            P(k,1) = flags.ID;
            P(k,2) = params.V0;
            P(k,3) = params.T0;
            P(k,4) = params.B0;
            P(k,5) = params.A0;
            P(k,6) = params.C0;
            P(k,7) = params.F0;
            P(k,8) = params.I0;
            P(k,9) = params.t0;
            P(k,10) = params.epsilon;
            P(k,11) = params.delta;
            P(k,12) = params.alphaL;
            P(k,13) = params.alphaV;
            P(k,14) = params.alphaT;
            P(k,15) = params.alphaB;
            P(k,16) = params.alphaC;
            P(k,17) = params.alphaI;
            P(k,18) = params.lambdaB;
            P(k,19) = params.lambdaC;
            P(k,20) = params.lambdaF;
            P(k,21) = params.lambdaI;
            P(k,22) = params.kappaF;
            P(k,23) = params.kappaI;
            P(k,24) = params.Ai/params.A0;
            P(k,25) = params.Fi/params.F0;
            P(k,26) = params.Ii/params.I0;
            
            T_scales(k,1) = flags.ID;
            T_scales(k,2) = params.delta*params.t0;
            T_scales(k,3) = params.t0;
            T_scales(k,4) = params.epsilon^(-1/3)*params.t0;
            T_scales(k,5) = params.epsilon^(-1)*params.t0;
            if params.alphaB<params.alphai
                T_scales(k,6) = -99;
            else
                T_scales(k,6) = -2/params.alphai/params.epsilon*log(params.epsilon)*params.t0;
            end
            T_scales(k,7) = tau(k);
        end
        
        T = array2table(P,'VariableNames',row_tits);
        writetable(T,'ndimparams.xlsx');
        
        T2 = array2table(T_scales,'VariableNames',t_tits);
        writetable(T2,'tscales.xlsx');
        
    case 1
        t_tits = {'ID','0','t1','t2','t3','t4','t5'};
        [C,S] = count_data;
        T = array2table(C,'VariableNames',t_tits);
        writetable(T,'datacount.xlsx');
    case 2
        row_tits = {'ID','tIstar','tImax','Imax','tBmax','Bmax','tAmax','Amax'};
         P = zeros(20,8);
         
        tau=[23;23;28;22;20;22;20;22;22;22;28;28;28;22;28;28;28;28;28;28];
       
        flags.paramset = 100;
        for k=1:20
            fprintf('Running data set %i\n',k);
            flags.ID = k;
            run_params;
            
            opt = get_optima([],params);
            
            P(k,1) = flags.ID;
            P(k,2) = (16*params.lambdaB/(params.lambdaI+params.lambdaB)^2)^(1/3);
            P(k,3) = params.t0*opt.tImax;
            P(k,4) = params.I0*opt.Imax;
            P(k,5) = params.t0*opt.tBmax;
            P(k,6) = params.B0*opt.Bmax;
            P(k,7) = params.t0*opt.tAmax;
            P(k,8) = params.A0*opt.Amax;
        end
        
        T = array2table(P,'VariableNames',row_tits);
        writetable(T,'optima.xlsx');
end
function [C,S,t,M] = count_data
flags.paramset = 100;
maxDATA = load('maxtimes.mat');
for k=1:20
    flags.ID = k;
    run_params;
    
    t(k,1) = params.delta*params.t0;
    t(k,2) = params.t0;
    t(k,3) = params.epsilon^(-1/3)*params.t0;
    t(k,4) = params.epsilon^(-1)*params.t0;
    if params.alphaB<params.alphai
        t(k,5) = -99;
    else
        t(k,5) = -2/params.alphai/params.epsilon*log(params.epsilon)*params.t0;
    end
    
    DATA = load(strcat(pwd,'\data\DATA_ID-',num2str(k),'.mat'));
    
    S(k,1) = k;
    C(k,1) = k;
    for j=1:5
        spot = find(DATA.X>t(k,j),1);
        if isempty(spot)
            S(k,j+1) = length(DATA.X);
        elseif spot==1
            S(k,j+1)=1;
        else
            S(k,j+1) = spot-1;
        end
        switch(j)
            case 1
                C(k,2) = sum(DATA.X<t(k,1));
            case 5
                if t(k,5) == -99
                    C(k,6) = sum(DATA.X>t(k,4));
                    C(k,7) = 0;
                else
                    C(k,6) = sum(DATA.X>t(k,4) & DATA.X<t(k,5));
                    C(k,7) = sum(DATA.X>=t(k,5));
                end
            otherwise
                C(k,j+1) = sum(DATA.X>t(k,j-1) & DATA.X<t(k,j));
        end
    end
    
    if k~= 8 || k~=10
        M(k,1) = k;
        M(k,2) = sum(DATA.X<=maxDATA.tImax(k));
        M(k,3) = sum(DATA.X>maxDATA.tImax(k) & DATA.X<=maxDATA.tBmax(k));
        M(k,4) = sum(DATA.X>maxDATA.tBmax(k) & DATA.X<=maxDATA.tAmax(k));
        M(k,5) = sum(DATA.X>maxDATA.tAmax(k));
    else
        M(k,1) = k;
        M(k,2) = sum(DATA.X<=maxDATA.tBmax(k));
        M(k,3) = sum(DATA.X>maxDATA.tBmax(k) & DATA.X<=maxDATA.tImax(k));
        M(k,4) = sum(DATA.X>maxDATA.tImax(k) & DATA.X<=maxDATA.tAmax(k));
        M(k,5) = sum(DATA.X>maxDATA.tAmax(k));
    end
    
end
end
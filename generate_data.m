function generate_data

T = readtable('Humoral_CytokineData_all2DoseData_logHumoral_logTFNy_CensoredIL_IncludesIL16_limitsIncluded_SutharInc');

for k=1:20
    ID_spots = find(T{20:end,2}==k)+19;
    censor = T{ID_spots,7};
    X = T{ID_spots,3};
    Y = T{ID_spots,4};
    
    X = X(X>=0);
    Y = Y(X>=0);
    censor = censor(X>=0);
    
    X = X(censor~=1);
    Y = Y(censor~=1);
    
    save(strcat(pwd,'\data\DATA_ID-',num2str(k),'.mat'),'X','Y');
end
end
if ~exist('flags','var')
    flags.paramset = 0;
end

%paramset
% -1 - artificial parameters for non-dimensional model
% 0 - population parameters from Chapin model
% 1 - Pfizer parameters
% 2 - Moderna parameters
% 100 - these are the different data sets from the literature given by
% patient IDs
switch flags.paramset
    case {0,-1}
        params.muLV = 9.1E-1;
        params.gammaL = 1.3E-4;
        params.gammaV = 7E-2;
        params.muTV = 4.98;
        params.gammaT = 5.5E-2;
        params.muTB = 9.8E-2;
        params.alphaBI = 2.4;
        params.sI = 1E3;
        params.gammaB = 7.1E-2;
        params.muBA = 0.51;
        params.gammaA = 4.2E-2;
        params.muCV = 2.2E-5;
        params.alphaCF = 1.4E-6;
        params.sF = 600;
        params.gammaC = 1E-2;
        params.muTF = 194.79;
        params.alphaFC = 1.3E-6;
        params.gammaF = 201.24;
        params.muTI = 2.07;
        params.alphaIB = 1.9E-3;
        params.gammaI = 2.7E-2;
        
        params.Ai = 22.2;
        params.Fi = 0;%5.3;
        params.Ii = 5.18;
        
    case 1
        params.muLV = 1.08;
        params.gammaL = 1.2E-4;
        params.gammaV = 7.34E-2;
        params.muTV = 5.96;
        params.gammaT = 5.81E-2;
        params.muTB = 9.84E-2;
        params.alphaBI = 2.47;
        params.sI = 1E3;
        params.gammaB = 0.128;
        params.muBA = 0.456;
        params.gammaA = 4.3E-2;
        params.muCV = 2.17E-5;
        params.alphaCF = 1.4E-6;
        params.sF = 600;
        params.gammaC = 1.04E-2;
        params.muTF = 193.63;
        params.alphaFC = 1.3E-6;
        params.gammaF = 201.42;
        params.muTI = 2.01;
        params.alphaIB = 2.2E-3;
        params.gammaI = 2.66E-2;
        
        params.Ai = 47.45;
        params.Fi = 0;
        params.Ii = 4.9;
    case 2
        params.muLV = 0.632;
        params.gammaL = 1.2E-4;
        params.gammaV = 6.4E-2;
        params.muTV = 10.53;
        params.gammaT = 5.16E-2;
        params.muTB = 9.52E-2;
        params.alphaBI = 2.34;
        params.sI = 1E3;
        params.gammaB = 5.08E-2;
        params.muBA = 0.736;
        params.gammaA = 4.16E-2;
        params.muCV = 2.26E-5;
        params.alphaCF = 1.4E-6;
        params.sF = 600;
        params.gammaC = 1.06E-2;
        params.muTF = 194.34;
        params.alphaFC = 1.3E-6;
        params.gammaF = 204.41;
        params.muTI = 14.3;
        params.alphaIB = 1.62E-3;
        params.gammaI = 2.76E-2;
        
        params.Ai = 292.53;
        params.Fi = 0;%5.3;
        params.Ii = 4.78;
    case 100
        T = readtable('paramdata');
        if ~isfield(flags,'ID')
            flags.ID = 1;
        end
        
        tabspot = find(T{:,1}==flags.ID);
        
        params.muLV = T{tabspot,4};
        params.gammaL = T{tabspot,5};
        params.gammaV = T{tabspot,6};
        params.muTV = T{tabspot,7};
        params.gammaT = T{tabspot,8};
        params.muTB = T{tabspot,9};
        params.alphaBI = T{tabspot,10};
        params.sI = T{tabspot,11};
        params.gammaB = T{tabspot,12};
        params.muBA = T{tabspot,13};
        params.gammaA = T{tabspot,14};
        params.muCV = T{tabspot,15};
        params.alphaCF = T{tabspot,16};
        params.sF = T{tabspot,17};
        params.gammaC = T{tabspot,18};
        params.muTF = T{tabspot,19};
        params.alphaFC = T{tabspot,20};
        params.gammaF = T{tabspot,21};
        params.muTI = T{tabspot,22};
        params.alphaIB = T{tabspot,23};
        params.gammaI = T{tabspot,24};
        params.Ai = T{tabspot,25};
        params.Ii = T{tabspot,26};
        params.Fi = 0;
        
end

%Asymptotic parameters
params.V0 = 1;
params.T0 = params.muTV/params.muLV*params.V0;
params.I0 = params.T0*params.muTI/params.muLV;
params.B0 = params.muTB*params.T0/params.muLV;
params.A0 = params.muBA*params.B0/params.gammaA;
params.C0 = params.muCV/params.muLV*params.V0;
params.F0 = params.muTF*params.T0/params.gammaF;

params.t0 = 1/params.muLV;

params.epsilon = params.gammaA/params.muLV;
params.delta = params.muLV/params.gammaF;

params.alphaL = params.gammaL/params.gammaA;
params.alphaV = params.gammaV/params.gammaA;
params.alphaB = params.gammaB/params.gammaA;
params.alphaT = params.gammaT/params.gammaA;
params.alphaI = params.gammaI/params.gammaA;
params.alphaC = params.gammaC/params.gammaA;

params.lambdaB = params.alphaBI*params.I0/params.sI/params.gammaA;
params.lambdaI = params.alphaIB*params.I0*params.muTB/params.muTI/params.gammaA;
params.lambdaC = params.alphaCF*params.muTF*params.muTV/params.sF/params.muLV/params.gammaA^2;
params.lambdaF = params.alphaFC*params.muCV/params.gammaA/params.muLV;

params.kappaI = params.I0/params.sI;
params.kappaF = params.F0/params.sF;

if flags.paramset == -1
    params.t0 = 1;
    params.V0 = 1;
    params.T0 = 1;
    params.B0 = 1;
    params.A0 = 1;
    params.C0 = 1;
    params.F0 = 1;
    params.I0 = 1;
    params.Ai = 1;
    params.Fi = 1;
    params.Ii = 1;
    
    params.epsilon = 0.01;
    params.delta = 0.01;
    
    params.alphaV = 2;%2;%1;
    params.alphaT = 1;%1;%2;
    params.alphaB = 3;%3;%3;
    params.alphaI = 1;
    params.alphaC = 1;
    
    params.lambdaB = 1;
    params.lambdaI = 1;
    
    params.alphaL = 0.01;
    params.lambdaC = 1;
    params.lambdaF = 1;
    params.kappaI = 0.01;
    params.kappaF = 0.01;
end

params.alphai = min(params.alphaT,params.alphaV);
params.alphaj = max(params.alphaT,params.alphaV);
params.Calphai = min(params.alphaV,params.alphaC);
params.Calphaj = max(params.alphaV,params.alphaC);
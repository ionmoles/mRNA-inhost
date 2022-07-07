function opt = get_optima(ts,params)

%Imax
% tol = 1E-6;
tmin = (16*params.lambdaB/(params.lambdaI+params.lambdaB)^2)^(1/3);
% tmin_old = tmin;

% done = 0;

% while ~done
%     Bmin = getBI(tmin,params);
%     tmin = (params.lambdaB+sqrt(params.lambdaB^2+2*params.lambdaI^3*Bmin^3*(params.lambdaI+params.lambdaB)))/(params.lambdaI*(params.lambdaI+params.lambdaB)*Bmin);
%     
%     if abs(tmin-tmin_old)<tol
%         done = 1;
%     else
%         tmin_old = tmin;
%     end
% end

tmaxI = fsolve(@(x) BI_fun(x,params),tmin);

opt.tImax = tmaxI*params.epsilon^(-1/3);
opt.Imax = (params.lambdaI+params.lambdaB)/2/params.lambdaB*tmaxI.^2-params.lambdaI/params.lambdaB*getBI(tmaxI,params);

%Bmax
alord = mink([params.alphaV,params.alphaT,params.alphaB],3);
ali = alord(1);
alj = alord(2);
alk = alord(3);

a = (alk-ali)*alj/(alk-alj)/ali;
b = (alj-ali)*alk/(alk-alj)/ali;

tmaxB = fsolve(@(x) 1-a*exp(-(alj-ali)*x)+b*exp(-(alk-ali)*x),1);
opt.tBmax = tmaxB/params.epsilon;
opt.Bmax = 1/params.epsilon^2*(1+params.lambdaB/params.lambdaI)/(params.alphaV-params.alphaT)*(exp(-params.alphaT*tmaxB)/(params.alphaB-params.alphaT)-exp(-params.alphaV*tmaxB)/(params.alphaB-params.alphaV)+(params.alphaV-params.alphaT)/((params.alphaB-params.alphaT)*(params.alphaB-params.alphaV))*exp(-params.alphaB*tmaxB));

%Amax
Aa(1) = 1/(params.alphaV-params.alphaT)/(params.alphaB-params.alphaT);
Aalpha(1) = params.alphaT;
Aa(2) = 1/(params.alphaT-params.alphaV)/(params.alphaB-params.alphaV);
Aalpha(2) = params.alphaV;
Aa(3) = 1/(params.alphaT-params.alphaB)/(params.alphaV-params.alphaB);
Aalpha(3) = params.alphaB;

tmaxA = fsolve(@(x) AB_fun(x,Aa,Aalpha),tmaxB);
opt.tAmax = tmaxA/params.epsilon;
opt.Amax = 0;
for k=1:3
    if Aalpha(k)==1
        opt.Amax = opt.Amax + (1+params.lambdaB/params.lambdaI)*Aa(k)*tmaxA*exp(-tmaxA);
    else
        opt.Amax = opt.Amax + (1+params.lambdaB/params.lambdaI)*Aa(k)/(1-Aalpha(k))*(exp(-Aalpha(k)*tmaxA)-exp(-tmaxA));
    end
end
opt.Amax = opt.Amax/params.epsilon^2;

end
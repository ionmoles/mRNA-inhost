function B=getBI(t,params)

bhat = (3*params.lambdaI+params.lambdaB)/3/(params.lambdaB+params.lambdaI);
Cu = gamma(bhat-1/3)/3/gamma(2/3);
chat = (params.lambdaB+params.lambdaI)/6;
zhat = chat*t.^3;
KM = hypergeom(bhat,4/3,zhat);
KU = kummerU(bhat,4/3,zhat);
KMp1 = hypergeom(bhat+1,4/3,zhat);
KUp1 = kummerU(bhat+1,4/3,zhat);

if t == 0
    B = 0;
else
    B = 1/params.lambdaI./t.*(1+3*bhat*((KMp1+Cu*(bhat-1/3)*KUp1)./(KM+Cu*KU)-1));
end
end
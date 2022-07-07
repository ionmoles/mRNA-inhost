function asdat = get_asymptotics(ts,flags,params,t_final)
%%%%t2 timescale

asdat.T0 = ts+exp(-ts)-1;
asdat.T1 = -(params.alphaV+params.alphaT)*(ts.^2/2-exp(-ts)-ts+1);

asdat.C0 = ts+exp(-ts)-1;
asdat.C1 = -(params.alphaV+params.alphaC)*(ts.^2/2-exp(-ts)-ts+1);

asdat.I0 = ts.^2/2-exp(-ts)-ts+1+params.Ii/params.I0;
aI = params.alphaV+params.alphaT+params.alphaI+params.lambdaI*(params.Ii/params.I0+2);
asdat.I1 = -(params.lambdaI/20*ts.^5-params.lambdaI/4*ts.^4+(params.lambdaI/3+aI/6)*ts.^3+(params.lambdaI*exp(-ts)-aI/2).*ts.^2-(params.lambdaI-params.alphaI*params.Ii/params.I0-aI)*ts-aI*(1-exp(-ts))+params.lambdaI/2*(1-exp(-2*ts)));

asdat.B0 = (ts.^2/2-exp(-ts)-ts+1);
a0 = params.alphaV+params.alphaT+params.alphaB-params.lambdaB*(params.Ii/params.I0+2);
asdat.B1 = params.lambdaB/20*ts.^5-params.lambdaB/4*ts.^4+(params.lambdaB/3-a0/6)*ts.^3+(a0/2+params.lambdaB*exp(-ts)).*ts.^2-(a0+params.lambdaB)*ts+a0*(1-exp(-ts))+params.lambdaB/2*(1-exp(-2*ts));

asdat.A0 = params.Ai/params.A0*ones(length(ts),1);
asdat.A1 = ts.^3/6-ts.^2/2+ts-1+exp(-ts)-params.Ai/params.A0*ts;

%%%%t3 timescale

that = params.epsilon^(1/3)*ts;

asdat.TI0 = that;
asdat.TI1 = -ones(length(that),1);

asdat.CI0 = that;
asdat.CI1 = -ones(length(that),1);

if flags.paramset<0
    kum_name = [pwd,'\Kummer\kummer-m',num2str(abs(flags.paramset)),'-',num2str(t_final),'-',num2str(length(ts)),'.mat'];
elseif flags.paramset==100
    kum_name = [pwd,'\Kummer\kummer-',num2str(abs(flags.paramset)),'-ID-',num2str(flags.ID),'-',num2str(t_final),'-',num2str(length(ts)),'.mat'];
else
    kum_name = [pwd,'\Kummer\kummer-',num2str(flags.paramset),'-',num2str(t_final),'-',num2str(length(ts)),'.mat'];
end

bhat = (3*params.lambdaI+params.lambdaB)/3/(params.lambdaB+params.lambdaI);
Cu = gamma(bhat-1/3)/3/gamma(2/3);
chat = (params.lambdaB+params.lambdaI)/6;
zhat = chat*that.^3;
if exist(kum_name,'file')==0
    fprintf('Kummer functions do not exist, computing...\n');
    KM = hypergeom(bhat,4/3,zhat);
    KU = kummerU(bhat,4/3,zhat);
    KMp1 = hypergeom(bhat+1,4/3,zhat);
    KUp1 = kummerU(bhat+1,4/3,zhat);
    save(kum_name,'KM','KU','KMp1','KUp1');
    fprintf('Kummer functions computed\n');
else
    fprintf('Kummer functions exist, loading...\n');
    temp = load(kum_name);
    KM = temp.KM;
    KU = temp.KU;
    KMp1 = temp.KMp1;
    KUp1 = temp.KUp1;
    fprintf('Kummer functions loaded\n');
end
asdat.BI0 = 1/params.lambdaI./that.*(1+3*bhat*((KMp1+Cu*(bhat-1/3)*KUp1)./(KM+Cu*KU)-1));
asdat.BI0(1) = 0;
Hhat = 2*params.lambdaI*cumtrapz(that,asdat.BI0)-chat*that.^3;
asdat.BI1 = -exp(-Hhat).*cumtrapz(that,(1+6*chat*that.*asdat.BI0).*exp(Hhat));
thatspot = find(abs(that-5)==min(abs(that-5)));
asdat.BI1(thatspot:end) = -(params.lambdaB+params.lambdaI)/params.lambdaI*that(thatspot:end);

asdat.II0 = (params.lambdaI+params.lambdaB)/2/params.lambdaB*that.^2-params.lambdaI/params.lambdaB*asdat.BI0;
asdat.II1 = -(params.lambdaI+params.lambdaB)/params.lambdaB*that-params.lambdaI/params.lambdaB*asdat.BI1;

asdat.AI0 = cumtrapz(that,asdat.BI0)+params.Ai/params.A0;
asdat.AI1 = cumtrapz(that,asdat.BI1);

%%%%t4 timescale
tau = params.epsilon*ts;

asdat.T0hat = (exp(-params.alphaT*tau)-exp(-params.alphaV*tau))/(params.alphaV-params.alphaT);
asdat.T1hat = (params.alphaT*exp(-params.alphaT*tau)-params.alphaV*exp(-params.alphaV*tau))/(params.alphaV-params.alphaT);

asdat.C0hat = (exp(-params.alphaC*tau)-exp(-params.alphaV*tau))/(params.alphaV-params.alphaC);
asdat.C1hat = (params.alphaC*exp(-params.alphaC*tau)-params.alphaV*exp(-params.alphaV*tau))/(params.alphaV-params.alphaC);

Chat = (params.alphaV-params.alphaT)/((params.alphaB-params.alphaT)*(params.alphaB-params.alphaV));
asdat.B0hat = (1+params.lambdaB/params.lambdaI)/(params.alphaV-params.alphaT)*(exp(-params.alphaT*tau)/(params.alphaB-params.alphaT)-exp(-params.alphaV*tau)/(params.alphaB-params.alphaV)+Chat*exp(-params.alphaB*tau));
asdat.B1hat = (1+params.lambdaB/params.lambdaI)/(params.alphaV-params.alphaT)*(params.alphaT*exp(-params.alphaT*tau)/(params.alphaB-params.alphaT)-params.alphaV*exp(-params.alphaV*tau)/(params.alphaB-params.alphaV)+params.alphaB*Chat*exp(-params.alphaB*tau));

asdat.I0hat = asdat.T0hat./(params.lambdaI*asdat.B0hat);
asdat.I1hat = asdat.T1hat./(params.lambdaI*asdat.B0hat) - asdat.T0hat.*asdat.B1hat./(params.lambdaI*asdat.B0hat.^2);

Iinf = (params.alphaB-params.alphai)/(params.lambdaI+params.lambdaB);

asdat.A0hat = exp(-tau).*cumtrapz(tau,asdat.B0hat.*exp(tau));%+params.epsilon^2*params.Ai/params.A0*exp(-tau);
asdat.A1hat = exp(-tau).*cumtrapz(tau,asdat.B1hat.*exp(tau));

%%%%t5 timescale
tauhat = tau+2/params.alphai*log(params.epsilon);

% Ttil = 1/(params.alphaV-params.alphaT)*(exp(-params.alphaT*tauhat)-params.epsilon^(2*(params.alphaV-params.alphaT)/params.alphaT)*exp(-params.alphaV*tauhat));
asdat.Ttil = 1/(params.alphaj-params.alphai)*exp(-params.alphai*tauhat);

asdat.Ctil = 1/(params.Calphaj-params.Calphai)*exp(-params.Calphai*tauhat);

%b = (1+params.lambdaB/params.lambdaI)/(params.alphaV-params.alphaT)*(1/(params.alphaB-params.alphaT)-params.epsilon^(2*(params.alphaV-params.alphaT)/params.alphaT)/(params.alphaB-params.alphaV)+Chat*params.epsilon^(2*(params.alphaB-params.alphaT)/params.alphaT));
% Btil = (1+params.lambdaB/params.lambdaI)/(params.alphaV-params.alphaT)*(exp(-params.alphaT*tauhat)/(params.alphaB-params.alphaT)-params.epsilon^(2*(params.alphaV-params.alphaT)/params.alphaT)*exp(-params.alphaV*tauhat)/(params.alphaB-params.alphaV)+Chat*params.epsilon^(2*(params.alphaB-params.alphaT)/params.alphaT)*exp(-params.alphaB*tauhat));
asdat.Btil = (1+params.lambdaB/params.lambdaI)/(params.alphaj-params.alphai)/(params.alphaB-params.alphai)*exp(-params.alphai*tauhat);

tauhatzerospot = find(abs(tauhat)==min(abs(tauhat)));
base_intB = trapz(tauhat(1:tauhatzerospot),asdat.Btil(1:tauhatzerospot));
G = params.alphaI*tauhat+params.lambdaI*(cumtrapz(tauhat,asdat.Btil)-base_intB);
base_intG = trapz(tauhat(1:tauhatzerospot),asdat.Ttil(1:tauhatzerospot).*exp(G(1:tauhatzerospot)));
asdat.Itil = exp(-G).*(cumtrapz(tauhat,asdat.Ttil.*exp(G))-base_intG)+Iinf*exp(-G);
end
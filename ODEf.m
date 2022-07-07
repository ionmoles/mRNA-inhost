function dy = ODEf(t,y,params,flags)

L = y(1);
V = y(2);
T = y(3);
B = y(4);
A = y(5);
C = y(6);
F = y(7);
I = y(8);

switch flags.modtype
    case 'full'

        dy(1) = -params.muLV*L - params.gammaL*L;
        dy(2) = params.muLV*L - params.gammaV*V;
        dy(3) = params.muTV*V - params.gammaT*T;
        dy(4) = params.muTB*T + params.alphaBI*(I/(I+params.sI))*B - params.gammaB*B;
        dy(5) = params.muBA*B - params.gammaA*A;
        dy(6) = params.muCV*V + params.alphaCF*(F/(F+params.sF))*C - params.gammaC*C;
        dy(7) = params.muTF*T - params.alphaFC*C*F - params.gammaF*F;
        dy(8) = params.muTI*T - params.alphaIB*I*B - params.gammaI*I;
    case 'scaled'
        dy(1) = -L - params.epsilon*params.alphaL*L;
        dy(2) = L - params.epsilon*params.alphaV*V;
        dy(3) = V - params.epsilon*params.alphaT*T;
        dy(4) = T + params.epsilon*params.lambdaB*I/(params.kappaI*I+1)*B - params.epsilon*params.alphaB*B;
        dy(5) = params.epsilon*(B - A);
        dy(6) = V + params.lambdaC*params.epsilon^2*params.delta*F/(params.kappaF*F+1)*C - params.epsilon*params.alphaC*C;
        dy(7) = 1/params.delta*(T - params.lambdaF*params.epsilon*params.delta*F*C - F);
        dy(8) = T - params.epsilon*params.lambdaI*I*B - params.epsilon*params.alphaI*I;
end

dy = dy(:);
end
function F= BI_fun(t,params)

B = getBI(t,params);

F = t - (params.lambdaB+sqrt(params.lambdaB^2+2*params.lambdaI^3*B^3*(params.lambdaI+params.lambdaB)))/(params.lambdaI*(params.lambdaI+params.lambdaB)*B);

end
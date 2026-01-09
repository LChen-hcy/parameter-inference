function S = Sfun_expression_LA_mean_fano(b)
%the negative log-likelihood function of telegraph model

global  pmdata  cellnumber meandata fanodata

kb_infimum = fanodata+meandata-1;
rho = exp(b(1))+kb_infimum;
sigma1 = meandata/rho*(rho+1-fanodata-meandata)/(fanodata-1); 
sigma0 = (rho-meandata)*sigma1/meandata;

N = length(pmdata);
pm = zeros(length(pmdata),1);
indices = find(pmdata);
%expression
for i = 1:1:length(indices)
    aa = sigma1+(indices(i)-1);
    bb = sigma1+sigma0+(indices(i)-1);
    x = -rho;
    y = 2*aa/(bb-x+sqrt((x-bb)^2+4*aa*x));
    r11 = y^2/aa+(1-y)^2/(bb-aa);
    pm(indices(i)) = (rho.^(indices(i)-1)./factorial(indices(i)-1)).*(gamma(sigma1+sigma0).*gamma(sigma1+(indices(i)-1)))./(gamma(sigma1).*gamma(sigma1+sigma0+(indices(i)-1))).*bb^(bb-1/2).*r11^(-1/2)*(y/aa)^aa*((1-y)/(bb-aa))^((bb-aa))*exp(x*y);
end

 pm_b = pm;
 pm_b(pm_b <= 0) =  1e-10;  

if min(pm_b) <= 0  ||  max(pm_b) >= 1.000000000 %remove calculation errors
    S = 1e+30;
else
        %the negative log-likelihood function
        S = 0;
        for i = 1:N
            S = S - cellnumber(i).*log(pm_b(i));
        end
end

end


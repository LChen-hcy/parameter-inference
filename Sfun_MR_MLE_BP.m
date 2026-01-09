function S = Sfun_MR_MLE_BP(b)
%the negative log-likelihood function of telegraph model

global  pmdata  cellnumber  meandata fanodata

kb_infimum = fanodata+meandata-1;
rho = exp(b(1))+kb_infimum;
sigma1 = meandata/rho*(rho+1-fanodata-meandata)/(fanodata-1);
sigma0 = (rho-meandata)*sigma1/meandata;
 lam2= 1;
N = length(pmdata);
pm = zeros(length(pmdata),1);
indices = find(pmdata);
%BP
for i = 1:1:length(indices)
    pm(indices(i)) = pBPi( indices(i)-1, sigma1, sigma0, rho, lam2);
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




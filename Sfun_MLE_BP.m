function S = Sfun_MLE_BP(b)
%the negative log-likelihood function of telegraph model

global  pmdata  cellnumber

rho =exp(b(3))+1; d = 1; sigma1 = exp(b(1)); sigma0 = exp(b(2)); lam2=1;

N = length(pmdata);
pm = zeros(length(pmdata),1);
indices = find(pmdata);
%BP
for i = 1:1:length(indices)
    pm(indices(i)) = pBPi(indices(i)-1, sigma1, sigma0, rho, lam2);
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




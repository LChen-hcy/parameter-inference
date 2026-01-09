function parfit_MLE_Exact
clc
close all
clear all

str={'estimating par','kon', ' koff', 'kb','kb/koff (bs)', 'time','AICc','HD'};

filename='result_MLE_Exact.xlsx';
data=xlsread('Embryonic_c57_distribution_1.xlsx');       %%Import data
xlswrite(filename,str,1,'A1')

for zushu=1:1:length(data(1,:))
    
    clearvars -except zushu data filename
    
    global  pmdata cellnumber
    
    N = 188;   % cell numbers
    
    pmdata = data(5:length(data(:,zushu)),zushu);
    pmdata = pmdata(~isnan(pmdata)); %remove NaN values
    
    S = length(pmdata);
    for i=1:1:S
        cellnumber(i)=N.*pmdata(i);
    end
    for i = 1:1:S
        cellnumber(i) = N.*pmdata(i);
    end
    
    tic
    % fit to telegraph model
    konini = 0.1 + 9.9*rand;
    kofini = 0.1 + 9.9*rand;
    kbini = 5 + 49*rand;
    ini = [log(konini)  log(kofini)  log(kbini)];
    [bminexpression, Sminexpression] = fminsearch(@Sfun_MLE_Exact,ini);  %call the negative log-likelihood function of telegraph model (FSP)
    kon_e = exp(bminexpression(1));
    koff_e = exp(bminexpression(2));
    kb_e = exp(bminexpression(3))+1;
    bs_e = kb_e/koff_e;
    toc
    time = toc;
    zushu
    
    k = 3;
    AICc_TM = 2*Sminexpression + 2*k + k*(k+1)/(N-k-1);
    
    rho =kb_e;  sigma1 = kon_e; sigma0 = koff_e;
    %expression
    for i = 1:1:length(pmdata)
        pm(i) = (rho.^(i-1)./factorial(i-1)).*(gamma(sigma1+sigma0).*gamma(sigma1+(i-1)))./(gamma(sigma1).*gamma(sigma1+sigma0+(i-1))).*hypergeom([sigma1+(i-1)],[sigma1+sigma0+(i-1)],-rho);
    end
    pm(pm <= 0) =  1e-10;

    err = 0;
    for i = 1:length(pmdata)
        err = (sqrt(pm(i))-sqrt(pmdata(i)))^2 + err;
    end
    HD = 1/sqrt(2)*sqrt(err);
    
    xlswrite(filename,[kon_e  koff_e  kb_e bs_e time AICc_TM HD],1,['B',num2str(zushu+1)])
    
    
end

end


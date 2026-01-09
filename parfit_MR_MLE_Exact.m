function parfit_MR_MLE_Exact
clc
close all
clear all

str={'estimating par','kon', ' koff', 'kb','kb/koff (bs)', 'time','AICc','HD'};

filename='result_expression_mean_fano.xlsx';
data=xlsread('Embryonic_c57_distribution_1.xlsx');       %%Import data
xlswrite(filename,str,1,'A1')

for zushu=1:1:length(data(1,:))
    
    clearvars -except zushu data filename
    
    global  pmdata cellnumber meandata fanodata
    
    N = 188;   % cell numbers
    
    pmdata = data(5:length(data(:,zushu)),zushu);
    pmdata = pmdata(~isnan(pmdata)); %remove NaN values
    
    meandata = 0;
    for i = 1:length(pmdata)
        meandata = meandata+(i-1).*pmdata(i);
    end
    u2data = 0;
    for i = 1:length(pmdata)
        u2data = u2data+(i-1)^2.*pmdata(i);
    end
    fanodata = u2data/meandata - meandata;
    
    S = length(pmdata);
    for i=1:1:S
        cellnumber(i)=N.*pmdata(i);
    end
    for i = 1:1:S
        cellnumber(i) = N.*pmdata(i);
    end
    
    tic 
    if fanodata <= 1
    kb_e = nan;
    kon_e = nan;
    koff_e = nan;
    bs_e = nan;
    Sminexpression_mean_fano = nan;
    else
    % fit to telegraph model
    kbini = 5 + 49*rand;
    ini = [log(kbini)];
    [bminexpression_mean_fano, Sminexpression_mean_fano] = fminsearch(@Sfun_MR_MLE_Exact,ini);  %call the negative log-likelihood function of telegraph model (FSP)
    kb_e = exp(bminexpression_mean_fano)+fanodata+meandata-1;
    kon_e = meandata/kb_e*(kb_e+1-fanodata-meandata)/(fanodata-1);
    koff_e = (kb_e-meandata)*kon_e/meandata;
    bs_e = kb_e/koff_e;
    end
    toc
    time = toc;
    zushu
    
    k = 1;
    AICc_TM = 2*Sminexpression_mean_fano + 2*k + k*(k+1)/(N-k-1);

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


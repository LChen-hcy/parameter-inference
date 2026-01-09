function parfit_MLE_Laplace
clc
close all
clear all

str={'estimating par','kon', ' koff', 'kb','kb/koff (bs)', 'time','AICc','HD'};

filename='result_MLE_Laplace.xlsx';
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
    [bminexpression_LA, Sminexpression_LA] = fminsearch(@Sfun_MLE_Laplace,ini);  %call the negative log-likelihood function of telegraph model (FSP)
    kon_e = exp(bminexpression_LA(1));
    koff_e = exp(bminexpression_LA(2));
    kb_e = exp(bminexpression_LA(3))+1;
    bs_e = kb_e/koff_e;
    toc
    time = toc;
    zushu
    
    k = 3;
    AICc_TM = 2*Sminexpression_LA + 2*k + k*(k+1)/(N-k-1);
    
    rho =kb_e;  sigma1 = kon_e; sigma0 = koff_e;
    %expression
    for i = 1:1:length(pmdata)
        aa = sigma1+(i-1);
        bb = sigma1+sigma0+(i-1);
        x = -rho;
        y = 2*aa/(bb-x+sqrt((x-bb)^2+4*aa*x));
        r11 = y^2/aa+(1-y)^2/(bb-aa);
        pm(i) = (rho.^(i-1)./factorial(i-1)).*(gamma(sigma1+sigma0).*gamma(sigma1+(i-1)))./(gamma(sigma1).*gamma(sigma1+sigma0+(i-1))).*bb^(bb-1/2).*r11^(-1/2)*(y/aa)^aa*((1-y)/(bb-aa))^((bb-aa))*exp(x*y);
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


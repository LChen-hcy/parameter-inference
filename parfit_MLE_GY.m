function parfit_MLE_GY
clc
close all
clear all

str={'estimating par','kon', ' koff', 'kb','kb/koff (bs)', 'time','AICc','HD'};

filename='result_MLE_GY.xlsx';
data=xlsread('Embryonic_c57_distribution_1.xlsx');       %%Import data
xlswrite(filename,str,1,'A1')

for zushu = 1:1:length(data(1,:))
    
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
    [bminBP, SminBP] = fminsearch(@Sfun_MLE_BP,ini);  %call the negative log-likelihood function of telegraph model 
    kon_e = exp(bminBP(1));
    koff_e = exp(bminBP(2));
    kb_e = exp(bminBP(3))+1;
    bs_e = kb_e/koff_e;
    toc
    time = toc;
    zushu

    k = 3;
    AICc_TM = 2*SminBP + 2*k + k*(k+1)/(N-k-1);
    
    %BP
    for i = 1:1:length(pmdata)
    pm(i) = pBPi(i-1, kon_e, koff_e, kb_e, 1);
    end
    
    err = 0;
    for i = 1:length(pmdata)
    err = (sqrt(pm(i))-sqrt(pmdata(i)))^2 + err;
    end
    HD = 1/sqrt(2)*sqrt(err);
    
    xlswrite(filename,[kon_e  koff_e  kb_e bs_e time AICc_TM HD],1,['B',num2str(zushu+1)])
    

end

end


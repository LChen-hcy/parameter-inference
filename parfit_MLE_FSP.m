function parfit_MLE_FSP
clc
close all
clear all

str={'estimating par','kon', ' koff', 'kb','kb/koff (bs)', 'time','AICc','HD'};

filename='result_MLE_FSP.xlsx';
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
    [bminFSP, SminFSP] = fminsearch(@Sfun_MLE_FSP,ini);  
    kon_e = exp(bminFSP(1));
    koff_e = exp(bminFSP(2));
    kb_e = exp(bminFSP(3))+1;
    bs_e = kb_e/koff_e;
    toc
    time = toc;
    zushu
    
    k = 3;
    AICc_TM = 2*SminFSP + 2*k + k*(k+1)/(N-k-1);
    
    rho =kb_e; d = 1; sigma1 = kon_e; sigma0 = koff_e;
    k=length(pmdata);
    for i=1:1:length(pmdata)
        if pmdata(k)>0
            break
        else
            k=length(pmdata)-i;
        end
    end
    vin=k-1;
    N = 3*(vin+1);
    
    % telegraph model FSP
    num = 2*N;
    Q = zeros(num);
    one = ones(num,1);
    for i = 1:N-1
        Q(i+1,i) = i*d;
        Q(i+N+1,i+N) = i*d;
        Q(i+N,i+N+1) = rho;
    end
    for i = 1:N
        Q(i,i+N) = sigma1;
        Q(i+N,i) = sigma0;
    end
    temp = Q*one;
    for i = 1:num
        Q(i,i) = -temp(i);
    end
    Qmod = Q;
    for i = 1:num
        Qmod(i,num) = 1;
    end
    vec = zeros(1,num);
    vec(num) = 1;
    ssd = vec/Qmod;
    dist = ssd(1:N)+ssd(N+1:2*N);
    NN = min(length(pmdata),N);
    pm = dist([1:NN]);
    pm(pm <= 0) =  1e-10;
    
    err = 0;
    for i = 1:NN
        err = (sqrt(pm(i))-sqrt(pmdata(i)))^2 + err;
    end
    HD = 1/sqrt(2)*sqrt(err);
    
    xlswrite(filename,[kon_e  koff_e  kb_e bs_e time AICc_TM HD],1,['B',num2str(zushu+1)])
    
end

end


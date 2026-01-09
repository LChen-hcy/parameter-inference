function S = Sfun_MR_MLE_FSP(b)
%the negative log-likelihood function of telegraph model

global  pmdata  cellnumber meandata fanodata

kb_infimum = fanodata+meandata-1;
rho = exp(b(1))+kb_infimum;
sigma1 = meandata/rho*(rho+1-fanodata-meandata)/(fanodata-1); 
sigma0 = (rho-meandata)*sigma1/meandata;
d = 1; 

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

 pm_b = pm;
 pm_b(pm_b <= 0) =  1e-10;  

if min(pm_b) <= 0  ||  max(pm_b) >= 1.000000000 %remove calculation errors
    S = 1e+30;
else
        %the negative log-likelihood function
        S = 0;
        for i = 1:NN
            S = S - cellnumber(i).*log(pm_b(i));
        end
end

end


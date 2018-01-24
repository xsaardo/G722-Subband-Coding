function rH = highBandADPCMDecoder(IH, step)

%% Constants
load('higherBandADPCMTable.mat');
N = length(IH);
delta_min = step/2;
beta = 127/128;

%% Convert codeword to table index
mH = zeros(length(IH),1);
mH(mod(IH,2) == 1) = 1;
mH(mod(IH,2) == 0) = 2;

%% initializations for HIGH
delta_H = ones(N,1);
dH = zeros(N,1);
div_H = zeros(N,1);
sHz = zeros(N,1);
sH = zeros(N,1);
rH = zeros(N,1);
pH = zeros(N,1);
aH1 = zeros(N,1); aH2 = zeros(N,1);
bH1 = zeros(N,1); bH2 = zeros(N,1); bH3 = zeros(N,1);
bH4 = zeros(N,1); bH5 = zeros(N,1); bH6 = zeros(N,1);

%%
for n = 7:length(IH)
    n
     %% quantizer adaptation
    div_H(n) = div_H(n - 1)*beta + WH(mH(n-1));
    
    if (div_H(n) < 0)
        div_H(n) = 0;
    end
    
    if (div_H(n) > 11)
        div_H(n) = 11;
    end
    
    delta_H(n) = 2^(div_H(n))*delta_min;
    
        %% 4 level inverse adaptive quantizer
    
    if (IH(n) == 0) || (IH(n) == 1)
        dH(n) = -Q2(mH(n))*delta_H(n);
    else
        dH(n) = Q2(mH(n))*delta_H(n);
    end
      %%
    sHp = aH1(n-1)*rH(n-1) + aH2(n-1)*rH(n-2);
    
    sHz(n)= bH1(n-1)*dH(n-1) + bH2(n-1)*dH(n-2) + bH3(n-1)*dH(n-3) ...
        + bH4(n-1)*dH(n-4) + bH5(n-1)*dH(n-5) + bH6(n-1)*dH(n-6);
    
    sH(n) = sHp + sHz(n);
    
        %% adaptive prediction
    rH(n) = sH(n) + dH(n);
    
    pH(n) = dH(n) + sHz(n);
    
    if abs(aH1(n-1))<=0.5
        f=4*(aH1(n-1));
    else f=2*sgn2(aH1(n-1));
    end
       
    pA = sgn2(pH(n))*sgn2(pH(n-1));
    pB = sgn2(pH(n))*sgn2(pH(n-2));
    
    aH1(n) = (1-2^-8)*aH1(n-1) + 3*2^-8*pA;
    aH2(n) = (1-2^-7)*aH2(n-1) + 2^-7*pB - 2^-7*f*pA;
    
    if aH2(n)>0.75
        aH2(n)=0.75;
    elseif aH2(n)<-0.75
        aH2(n)=-0.75;
    end
    
    if aH1(n)>(1-2^-4 - aH2(n))
        aH1(n)=1-2^-4-aH2(n);
    elseif aH1(n)<-(1-2^-4 - aH2(n))
        aH1(n)=-(1-2^-4 -aH2(n));
    end
    
    bH1(n) = (1-2^(-8))*bH1(n-1) + 2^(-7)*sign(dH(n))*sgn2(dH(n-1));
    bH2(n) = (1-2^(-8))*bH2(n-1) + 2^(-7)*sign(dH(n))*sgn2(dH(n-2));
    bH3(n) = (1-2^(-8))*bH3(n-1) + 2^(-7)*sign(dH(n))*sgn2(dH(n-3));
    bH4(n) = (1-2^(-8))*bH4(n-1) + 2^(-7)*sign(dH(n))*sgn2(dH(n-4));
    bH5(n) = (1-2^(-8))*bH5(n-1) + 2^(-7)*sign(dH(n))*sgn2(dH(n-5));
    bH6(n) = (1-2^(-8))*bH6(n-1) + 2^(-7)*sign(dH(n))*sgn2(dH(n-6));
    
    %% 
    rH(n) = dH(n) + sH(n);
    
end

end
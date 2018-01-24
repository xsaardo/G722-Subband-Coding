function rL = lowBandADPCMDecoder(IL,step,mode)

%% Constants
load('lowerBandADPCMTable.mat');
N = length(IL);
beta = 127/128;
delta_min = step/2;

QL6 = lowerBandADPCMTable.QL6;
QL5 = lowerBandADPCMTable.QL5;
QL4 = lowerBandADPCMTable.QL4;
WL = lowerBandADPCMTable.WL;
pos_code = lowerBandADPCMTable.PositiveCodes;
neg_code = lowerBandADPCMTable.NegativeCodes;

%% Convert codeword to table index
mL = zeros(length(IL),1);
for i = 1:length(IL)
    [a,~] = find(pos_code == IL(i));
    [b,~] = find(neg_code == IL(i));
    if (isempty(a))
        mL(i) = b;
    else
        mL(i) = a;
    end
end

%% initializations
delta_L = ones(N,1);
dLt = zeros(N,1);
div_L = zeros(N,1);
sLz = zeros(N,1);
sL = zeros(N,1);
rLt = zeros(N,1);
pLt = zeros(N,1);
aL1 = zeros(N,1); aL2 = zeros(N,1);
bL1 = zeros(N,1); bL2 = zeros(N,1); bL3 = zeros(N,1);
bL4 = zeros(N,1); bL5 = zeros(N,1); bL6 = zeros(N,1);
dL = zeros(N,1);
rL = zeros(N,1);

%%
for n = 7:length(IL)
    n
    %% quantizer adaptation
    div_L(n) = div_L(n - 1)*beta + WL(mL(n-1));
    
    if (div_L(n) < 0)
        div_L(n) = 0;
    end
    
    if (div_L(n) > 9)
        div_L(n) = 9;
    end
    
    delta_L(n) = 2^(div_L(n)+2)*delta_min;
    
    %% 60 and 15 level inverse adaptive quantizer
    if((IL(n) >= 4 && IL(n) <= 31) || (IL(n) == 62) || (IL(n) == 63))
        dLt(n) = -QL4(mL(n))*delta_L(n);
        % Mode selection
        if (mode == 1)
            dL(n) = -QL6(mL(n))*delta_L(n);
        elseif (mode == 2)
            dL(n) = -QL5(mL(n))*delta_L(n);
        else
            dL(n) = -QL4(mL(n))*delta_L(n);
        end
    else
        dLt(n) = QL4(mL(n))*delta_L(n);
        % Mode selection
        if (mode == 1)
            dL(n) = QL6(mL(n))*delta_L(n);
        elseif (mode == 2)
            dL(n) = QL5(mL(n))*delta_L(n);
        else
            dL(n) = QL4(mL(n))*delta_L(n);
        end
    end
    
    %%
    sLp = aL1(n-1)*rLt(n-1) + aL2(n-1)*rLt(n-2);
    
    sLz(n)= bL1(n-1)*dLt(n-1) + bL2(n-1)*dLt(n-2) + bL3(n-1)*dLt(n-3) ...
        + bL4(n-1)*dLt(n-4) + bL5(n-1)*dLt(n-5) + bL6(n-1)*dLt(n-6);
    
    sL(n) = sLp + sLz(n);
    
    %% adaptive prediction
    rLt(n) = sL(n) + dLt(n);
    
    pLt(n) = dLt(n) + sLz(n);
    
    if abs(aL1(n-1))<=0.5
        f=4*(aL1(n-1));
    else f=2*sgn2(aL1(n-1));
    end
    
    pA = sgn2(pLt(n))*sgn2(pLt(n-1));
    pB = sgn2(pLt(n))*sgn2(pLt(n-2));
    
    aL1(n) = (1-2^-8)*aL1(n-1) + 3*2^-8*pA;
    aL2(n) = (1-2^-7)*aL2(n-1) + 2^-7*pB - 2^-7*f*pA;
    
    if aL2(n)>0.75
        aL2(n)=0.75;
    elseif aL2(n)<-0.75
        aL2(n)=-0.75;
    end
    
    if aL1(n)>(1-2^-4 - aL2(n))
        aL1(n)=1-2^-4-aL2(n);
    elseif aL1(n)<-(1-2^-4 - aL2(n))
        aL1(n)=-(1-2^-4 -aL2(n));
    end
    
    bL1(n) = (1-2^(-8))*bL1(n-1) + 2^(-7)*sign(dLt(n))*sgn2(dLt(n-1));
    bL2(n) = (1-2^(-8))*bL2(n-1) + 2^(-7)*sign(dLt(n))*sgn2(dLt(n-2));
    bL3(n) = (1-2^(-8))*bL3(n-1) + 2^(-7)*sign(dLt(n))*sgn2(dLt(n-3));
    bL4(n) = (1-2^(-8))*bL4(n-1) + 2^(-7)*sign(dLt(n))*sgn2(dLt(n-4));
    bL5(n) = (1-2^(-8))*bL5(n-1) + 2^(-7)*sign(dLt(n))*sgn2(dLt(n-5));
    bL6(n) = (1-2^(-8))*bL6(n-1) + 2^(-7)*sign(dLt(n))*sgn2(dLt(n-6));
    
    %%
    rL(n) = dL(n) + sL(n);
    
end

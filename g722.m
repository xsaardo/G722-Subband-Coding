% input speech signal
[x,Fs] = audioread('man1_wb.wav');

% resample to 16KHz
x = resample(x,16000/Fs,1);
%x = x(1:16000);

% Quantize to 14 bits
low = min(x);
high = max(x);
step = (high - low)/(2^14 - 1);
quantization_steps = low:step:high;
x = discretize(x,quantization_steps);

% Normalize
x = x/max(x);

% coefficients for h0-h11, h23-h12
qmf = [0.366211E-03, -0.134277E-02, -0.134277E-02, 0.646973E-02, 0.146484E-02, -0.190430E-01, 0.390625E-02, 0.441895E-01, -0.256348E-01, -0.982666E-01 0.116089E+00, 0.473145E+00];
h0 = [qmf fliplr(qmf)];

%% QMF filtering
i = 1:length(h0);
h1 = (-1).^i.*h0;

xL = downsample(conv(x,h0),2);
xH = downsample(conv(x,h1),2);

% subplot(2,1,1);
% plot(xL);
% subplot(2,1,2);
% plot(xH);

%%
N = length(xL);

decision_levels = [adapQuantTable(:,1);inf];
pos_code = adapQuantTable(:,6);
neg_code = adapQuantTable(:,5);
QL4 = adapQuantTable(:,3);
WL = adapQuantTable(:,4);

%% initializations
eL = zeros(N,1);
delta_L = ones(N,1);
dLt = zeros(N,1);
div_L = zeros(N,1);
Ilt = zeros(N,1)*15;
IL = zeros(N,1)*63;
code = ones(N,1);
sLz = zeros(N,1);
sL = zeros(N,1);
rLt = zeros(N,1);
pLt = zeros(N,1);
aL1 = zeros(N,1); aL2 = zeros(N,1);
bL1 = zeros(N,1); bL2 = zeros(N,1); bL3 = zeros(N,1);
bL4 = zeros(N,1); bL5 = zeros(N,1); bL6 = zeros(N,1);

%% constants
beta = 127/128;
delta_min = step/2;

%%
for n = 7:length(xL)
    n
     %% quantization adaption
    div_L(n) = div_L(n - 1)*beta + WL(code(n-1));
    
    if (div_L(n) < 0)
        div_L(n) = 0;
    end
    
    if (div_L(n) > 9)
        div_L(n) = 9;
    end
    
    delta_L(n) = 2^(div_L(n)+2)*delta_min;
    
    %%
    sLp = aL1(n-1)*rLt(n-1) + aL2(n-1)*rLt(n-2);
    
    sLz(n)= bL1(n-1)*dLt(n-1) + bL2(n-1)*dLt(n-2) + bL3(n-1)*dLt(n-3) ...
        + bL4(n-1)*dLt(n-4) + bL5(n-1)*dLt(n-5) + bL6(n-1)*dLt(n-6);
    
    sL(n) = sLp + sLz(n);
    
    %% 60 level quantization, truncate last 2 bits
    eL(n) = xL(n)  - sL(n);
    if (eL(n) < 0)
        code(n) = discretize(-eL(n), decision_levels*delta_L(n));
        IL(n) = neg_code(code(n));
        dLt(n) = -QL4(code(n))*delta_L(n);
    else
        code(n) = discretize(eL(n), decision_levels*delta_L(n));
        IL(n) = pos_code(code(n));
        dLt(n) = QL4(code(n))*delta_L(n);
    end
    
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
   
end

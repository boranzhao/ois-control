%*************************************************************************
%  
%  Hinf and mu-synthesis controller design for miniaturized optical image stabilizers 
%             ----------------------------------------
%
% By Pan Zhao, Control Engineering Lab, UBC, BC, Canada
% Under supervision of Dr. Ryozo Nagamune.
% Revision: Feb 25, 2018
% Creation: Jan 21, 2016.

% Function: Main function for designing a Hinf and mu-synthesis
% Reference: Zhao, R. Nagamune and M. Chiao, \Multiple parameter-dependent robust control of miniaturized optical
% image stabilizers," accepted by Control Engineering Practice, 2018.

clear;close all;
%% Plot Frequency Response 

% Load filter data
load('Filters.mat')
%% parameters 
Ktype = 2;              %  1 for H infinity control; 2 for mu control
PLT_INDEX = 3;          %1-5 for plant 1-5
FREQ_NQST = PLT_INDEX;  %0 for model; 1-5 for data 1-5;
MODE1_UNCERT = 1;       % 1: consider uncertainty brought by w2n for mode 1; 0: ignore uncertainty brought by w2n for mode 1
BSFILTER = 1;           %1: include band stop filter; 0: without band stop filter 
BSFILTER_MOD1 = 0;      % whether to include the notch filter for Mode 1
Ts = 5e-4;              % sample time for model discretization 
fig1 = figure;
%% Load and Plot frequency response data for all the plants
% bode & nyquist plot setting
P = bodeoptions; % Set phase visiblity to off and frequency units to Hz in options
P.FreqUnits = 'Hz'; % Create plot with the options specified by P
    
Pn = nyquistoptions;
Pn.FreqUnits = 'Hz'; 

ReadFreqResp();

subplot(2,1,1)
goodplot;
subplot(2,1,2)
goodplot([7 7]);
% print -painters -dpdf -r150 MeasFR.pdf

%% Estimate a model based on measured frequency response data. 
DCgain = 10^(12.29/20); %12.41;
% zeta1 = 5e-3; % for mu controller
w10 = 63.84*2*pi;
w20 = 101.7*2*pi; %104.8
% 2nd-order model for Mode 1
if PLT_INDEX == 1
    w20 = 97.11*2*pi; %104.8
    w10 = 62.35*2*pi;
    K1 = -DCgain/10;
    zeta1 = 5e-3
elseif PLT_INDEX == 2
    w20 = 101.6*2*pi; %104.8
    w10 = 62.84*2*pi; %w20/101.7*63.84
    K1 = DCgain/20;
elseif PLT_INDEX == 3
    w20 = 101.7*2*pi; %104.8
    w10 = 63.84*2*pi;
    zeta1 = 5e-3;
    K1 = DCgain/10;%20 successful for nomial plant, 50
elseif PLT_INDEX == 4
    w20 = 104.9*2*pi; %104.8
    w10 = 67.33*2*pi;
    K1 = DCgain/20;%20 successful for nomial plant, 50
    zeta1 = 10e-3;
elseif PLT_INDEX == 5  
    w20 = 106.2*2*pi; %104.8
    w10 = 69.5*2*pi;
    K1 = DCgain/100;
    zeta1 = 10e-3;
end
for iii = 1
% for 4th-order model for Mode 1
% if PLT_INDEX == 1
%     w20 = 97.11*2*pi; %104.8
%     w10 = 62.35*2*pi;
%     K1 = -DCgain/10;
%     zeta1 = 5e-3
% elseif PLT_INDEX == 2
%     w20 = 101.6*2*pi; %104.8
%     w10 = 62.84*2*pi; %w20/101.7*63.84
%     K1 = DCgain/20;
% elseif PLT_INDEX == 3
%     w20 = 101.7*2*pi; %104.8
%     w10 = 63.84*2*pi;
%     zeta1 = 1.5e-2;
%     K1 = DCgain/8;%20 successful for nomial plant, 50
% elseif PLT_INDEX == 4
%     w20 = 104.9*2*pi; %104.8
%     w10 = 67.33*2*pi;
%     K1 = DCgain/20;%20 successful for nomial plant, 50
%     zeta1 = 10e-3;
% elseif PLT_INDEX == 5  
%     w20 = 106.2*2*pi; %104.8
%     w10 = 69.5*2*pi;
%     K1 = DCgain/100;
%     zeta1 = 10e-3;
% end
end

zeta1 = 1e-3; % for LPV controller
switch Ktype 
    case 1
        wn = w20*1.00;
    case 2
%         wn = ureal('wn',w20,'perc',5,'AutoSimplify','full'); % for mu
%         controller
%         wn = ureal('wn',w20,'range',[96.615  106.785]*2*pi,'AutoSimplify','full');
        wn = ureal('wn',w20,'range',[95 108]*2*pi,'AutoSimplify','full'); % for the LPV model
%         wn = ureal('wn',w20,'PlusMinus',[96.615  106.785]*2*pi-w20,'AutoSimplify','full');
%         K1 = ureal('K1',0,'PlusMinus',0.15,'AutoSimplify','full'); % for
%         mu controller
        K1 = ureal('K1',0,'PlusMinus',0.2,'AutoSimplify','full'); % for LPV controller
end  
% Mode 1
% for i = 1
% Model 1: the first mode
% wn1n= w10*0.99;
% zeta1n = 1e-2;
% 
% wn1d_nom = w10;
% wn1d = wn1d_nom;
% K1 = wn1n^2/wn1d^2;
% zeta1d = 1e-3;
% Model1 = tf([1 2*zeta1n*wn1n wn1n^2],[1 2*zeta1d*wn1d wn1d^2]);%*tf(wn1d^2,[1 2*zeta1d*wn1d wn1d^2])
% end
%Mode 1: the first mode

%%4th-order
% K1 = DCgain/500;%50
% zeta1 = 2.23e-3; %1.5e-2;
% zeta1 = 5e-3; %1.5e-2;
if MODE1_UNCERT == 1
%     wn1 = w10*wn/w20; 
%      delta = ureal('delta',0,'PlusMinus',[-3 3],'AutoSimplify','full'); %
%      for mu controller
     wn1 = ureal('wn1',0,'range',[61 71]*2*pi,'AutoSimplify','full') ; % for LPV controller
else
    wn1 = w10;
end
Model1 = K1*tf(wn1^2,[1 2*zeta1*wn1 wn1^2]);
% Mode 2  
K2 = DCgain; %% (DCgain-K1) This may cause significant change of controller dynamics ; %
% K2 = DCgain;
%  K2 = DCgain/K1;
%K_nom = w20^2*K2;
zeta2 = 5e-4;
% zeta2 = 1.5e-3;
Model2 = K2*(w20^2/wn^2)*tf(wn^2,[1 2*zeta2*wn wn^2]); 
% Mode 3&4
for i=1
% % Model 3 
% wn3n = 243.4*2*pi;
% zeta3n = 2e-3;
% 
% wn3d = 254.4*2*pi;
% K3 = wn3d^2/wn3n^2;
% zeta3d = 2e-3;
% Model3 = tf(K3*[1 2*zeta3n*wn3n wn3n^2],[1 2*zeta3d*wn3d wn3d^2]);
% 
% % Model 4
% wn4n = 257.4*2*pi;
% zeta4n = 2e-3;
% wn4d = 271.3*2*pi;
% % wn4d = 263.3*2*pi;
% K4 = wn4d^2/wn4n^2;
% zeta4d = 2e-3;
% Model4 = tf(K4*[1 2*zeta4n*wn4n wn4n^2],[1 2*zeta4d*wn4d wn4d^2]);
end

% time delay estimation
DelayT = 21.86/(360*52.87);
[numD,denD] = pade(DelayT,1);
DelayMod = tf(numD,denD);

% final model 
Model = (Model1+Model2)*DelayMod; %Model1 Model2+
if Ktype == 1
    sys = ss(Model); 
else
    sys = Model;
    Model = sys.NominalValue;
end
ModelSS = ss(Model);
%RedOrder = order(sys);
%% frequency response of the model (Hinf) or nominal model( mu)
freqM = (1:0.1:400)';
[magM,phaseM] = bode(Model,freqM*2*pi);
phaseM = phaseM(:);
% phase4 = phase4-freq*2*pi*DelayT*180/pi;
for i= 1:size(phaseM,1)
    if phaseM(i,1)>360
        phaseM(i,1) = phaseM(i,1)-720;
    elseif phaseM(i,1)>40
        phaseM(i,1) = phaseM(i,1)-360;
    end    
end
magM = 20*log10(magM(:));
figure(fig1)
subplot(2,1,1);semilogx(freqM,magM,'-.m','Linewidth',2);
xlabel('Frequency (Hz)','Fontsize',13);ylabel('Magnitued (dB)','Fontsize',13);
legend('Exp 1','Exp 2','Exp 3','Exp4','Exp5','Model')
grid on;hold off;set(gca,'Fontsize',12);
subplot(2,1,2);semilogx(freqM,phaseM,'-.m','Linewidth',2);
xlabel('Frequency (Hz)','Fontsize',13);ylabel('Phase (deg)','Fontsize',13);
legend('Exp 1','Exp 2','Exp 3','Exp4','Exp5','Model');
set(gca,'Fontsize',12);grid on;
title('Frequency response data and of nomial model')
% semilogx(freq,-180*ones(size(phase)),'k')
%% frequency response of randomly sampled models from the uncertain model
% if Ktype == 2
% fig2 = figure(); 
% systmp = usample(sys,15);
%     for kk = 1:10
%         [magM,phaseM] = bode(systmp(:,:,kk,1),freqM*2*pi);
%         magM = 20*log10(magM(:));
%         phaseM = phaseM(:);
%         for i= 1:size(phaseM,1)
%             if phaseM(i,1)>370
%                 phaseM(i,1) = phaseM(i,1)-720;
%             elseif phaseM(i,1)>10
%                 phaseM(i,1) = phaseM(i,1)-360;
%             end
%         end
%         subplot(2,1,1);semilogx(freqM,magM,'Linewidth',2);hold on;grid on;
%         ylabel('Magnitude (dB)')
%         subplot(2,1,2);semilogx(freqM,phaseM,'Linewidth',2);hold on;grid on;
%         xlabel('Frequency (Hz)')
%         ylabel('Phase (deg)')
%     end
% end
% goodplot;
% print -painters -dpdf -r150 MeasFR.pdf

%% frequency response of sampled models with specific parameters from the uncertain model
fig2 = figure(); 
% [systmp, sampleValues] = usample(sys,8);
color = {'-.g','b--','r','c--','-.k','m--','r--','b:'};
lineWidth = [1.5 1 1 1 1 1 1 1];
markerSize = [3 1 2 3 2 2 2 2];
wnVec = [95 97 100 105 108]*2*pi;
wn1Vec =[67 63 65 67 71]*2*pi;
K1Vec = [-0.1 -0.2 0.2 0.1 0];

legendTxt = cell(1,5);
for kk = 1:5
    % get the sampled values
%     wn1 = sampleValues(kk).wn1;
%     wn2 = sampleValues(kk).wn;
%     K1 = sampleValues(kk).K1;
%     sysTest = systmp(:,:,kk,1);

    wn1 = wn1Vec(kk);
    wn2 = wnVec(kk);
    K1 = K1Vec(kk);
    sysTest = usubs(sys,'wn1',wn1,'wn',wn2,'K1',K1);
    
    [magM,phaseM] = bode(sysTest,freqM*2*pi);
    magM = 20*log10(magM(:));
    phaseM = phaseM(:);
%     legendTxt{kk} = sprintf('$(\\omega_1,\\omega_2) = (%.1f, %.1f)$ Hz, $K_1= %.2f$',wn1/2/pi, sampleValues(kk).wn/2/pi, sampleValues(kk).K1);
    legendTxt{kk} = sprintf('$(\\omega_1,\\omega_2,K1) = (%.0f, %.0f, %.1f)$',wn1/2/pi, wn2/2/pi, K1);
  
    for i= 1:size(phaseM,1)
        if phaseM(i,1)>370
            phaseM(i,1) = phaseM(i,1)-720;
        elseif phaseM(i,1)>20
            phaseM(i,1) = phaseM(i,1)-360;
        end
    end
    subplot(2,1,1);semilogx(freqM,magM,color{kk},'Linewidth',lineWidth(kk),'MarkerSize',markerSize(kk));hold on;grid on;
    ylabel('Magnitude (dB)')
    subplot(2,1,2);semilogx(freqM,phaseM,color{kk},'Linewidth',lineWidth(kk),'MarkerSize',markerSize(kk));hold on;grid on;
    xlabel('Frequency (Hz)')
    ylabel('Phase (deg)')
end
    
subplot(2,1,1)
xlim([1e1 1e3])
ylim([-50 100])
legend(legendTxt,'interpreter','latex');
goodplot;
subplot(2,1,2)
xlim([1e1 1e3])
ylim([-400 50])
legend(legendTxt,'interpreter','latex');
goodplot([7 7]);
% print -painters -dpdf -r150 sampleModelFR.pdf

return;
%% weighting function
Weinv = makeweight(0.028,750,3);%(0.001,20,2);
WTinv = makeweight(3,100,0.001);
We = 1/5*ss(inv(Weinv));
WT = inv(WTinv);
%penalize the input around the 2nd mode
% Wu = 1/100*1/5*ss(11.5*tf([1/(w20*0.5) 1],[1/(w20*0.7) 1])*tf([1/(w20*1.4) 1],[1/(w20*1) 1])); 
Wu = 0.023;

if Ktype == 1 % H infinity controller design
    for i1 =1
%         Wu = 10;     
        n = size(sys.A,1);
        nws = 1;
        nwu = 2;  
        if order(ss(Wu)) == 0
            sysG = ss([sys.A zeros(size(sys.A,1),1); -We.B*sys.C We.A], [zeros(size(sys.A,1),1) sys.B; We.B -We.B*sys.D],... % consider constant Wu
                     [-We.D*sys.C We.C;zeros(1,size(sys.A,1)) 0; -sys.C 0], [We.D -We.D*sys.D+1e-3;0 Wu; 1 -sys.D+1e-10]);
        elseif order(Wu) > 0
            alp = 0;
            sysGA = [sys.a zeros(n,nws) zeros(n,nwu);
                 -We.b*sys.c We.A zeros(nws,nwu);
                 zeros(nwu,n) zeros(nwu,nws) Wu.a]; sysGA = sysGA+alp*eye(size(sysGA));
             sysGB =  [zeros(n,1) sys.B; We.B -We.B*sys.D;zeros(nwu,1) Wu.B];
             sysGC = [-We.d*sys.c We.c zeros(1,nwu);
                 zeros(1,n) zeros(1,nws) Wu.c; 
                 -sys.C 0 zeros(1,nwu)];
             sysGD = [We.d -We.d*sys.d+1e-3;0 Wu.d;1 -sys.d];
            sysG = ss(sysGA,sysGB,sysGC,sysGD); %consider dynamic Wu
        end
        Wu = ss(Wu);
            % sysG = ss([sys.A zeros(size(sys.A,1),1); -We.B*sys.C We.A], [zeros(size(sys.A,1),1) sys.B; We.B -We.B*sys.D],...
        %[-We.D*sys.C We.C;-sys.C 0], [We.D -We.D*sys.D+1e-3;1 -sys.D]); % igore Wu
        RedOrder = order(sysG);
%         [sysGmin] = balreal(sysG);
        [Khinf1,clp,Gam] = hinfsyn(sysG,1,1,'method','lmi'); Gam
        break;
%         Khinf1.a = Khinf1.a-alp*eye(size(Khinf1.a));
        %reduce the order of the model
        KhinfR = balred(Khinf1,RedOrder);
        figure;
        sigma(Khinf1,KhinfR)
        legend('Full-order model','Reduced-order model');
        sysS = 1/(1+sys*KhinfR);
        %Plot frequency response
        figure
        bodeplot(Weinv,inv(Wu),sysS,KhinfR*sysS,P)%Khinf,1-sysS,,sys*sysK/(1+sys*Khinf),P
        legend('Weighting function','Weighting function u','Senstivity','Control Input')%,,'Complementary sensitivity','Controller','Closed-loop freq resp'
        grid on;                     
        %% Continuous-Discrete model transformation
        KhinfRD = c2d(KhinfR,Ts,'prewarp',90*2*pi)
        BSfilterD = c2d(BSfilter,Ts,'prewarp',252*2*pi)
        BSfilterMod1D = c2d(BSfilterMod1,Ts,'prewarp',65*2*pi)        
        Khinf = KhinfR*BSfilter; %
        % Determine whether to add the band stop filter 
        if BSFILTER == 1
            Khinf = KhinfR * BSfilter;
            KhinfD = KhinfRD * BSfilterD;
        elseif BSFILTER == 0
            Khinf = KhinfR;
            KhinfD = KhinfRD;
        end
        KhinfD = balreal(KhinfD);
        figure;
        bodeplot(Khinf,KhinfD,P);
        legend('Continous','Discrete');
        [numKhinfD,denKhinfD] = tfdata(KhinfD,'v');    
        [numKhinf,denKhinf] = tfdata(Khinf,'v'); 
        sysK = Khinf;
    end    
elseif Ktype == 2 % mu controller
    for foldthecode =1
%        Wu = ss([],[],[],5);
%        Wu = makeweight(0.1,500,6);
        mws = 1;
        n = size(sys.A,1);
        nws = 1;
        nwu = 2;
        if order(ss(Wu)) == 0
            sysG = ss([sys.A zeros(size(sys.A,1),1); -We.B*sys.C We.A], [zeros(size(sys.A,1),1) sys.B; We.B -We.B*sys.D],... % consider constant Wu
                     [-We.D*sys.C We.C;zeros(1,size(sys.A,1)) 0; -sys.C 0], [We.D -We.D*sys.D+1e-3;0 Wu; 1 -sys.D+1e-10]);
        elseif order(Wu) > 0
            sysGA = [sys.a zeros(n,nws) zeros(n,nwu);
                 -We.b*sys.c We.A zeros(nws,nwu);
                 zeros(nwu,n) zeros(nwu,nws) Wu.a];
            sysGB =  [zeros(n,1) sys.B; We.B -We.B*sys.D;zeros(nwu,1) Wu.B];
            sysGC = [-We.d*sys.c We.c zeros(1,nwu);
                 zeros(1,n) zeros(1,nws) Wu.c; 
                 -sys.C 0 zeros(1,nwu)];
            sysGD = [We.d -We.d*sys.d+1e-3;0 Wu.d;1 -sys.d];
            sysG = ss(sysGA,sysGB,sysGC,sysGD); % consider dynamic Wu
        end
        %sysG = ss([sys.A zeros(size(sys.A,1),1); -We.B*sys.C We.A], [zeros(size(sys.A,1),1) sys.B; We.B -We.B*sys.D],...
            %[-We.D*sys.C We.C;-sys.C 0], [We.D -We.D*sys.D+1e-3;1 -sys.D]); % igore Wu
        opt = dkitopt('FrequencyVector',logspace(-2,2.6,250)*2*pi,'NumberOfAutoIterations',20,'MixedMU','on','AutoScalingOrder',[5 2]); %
        [Kmu1,clp,bnd,DKinfo] = dksyn(sysG,1,1,opt);       
%         reduce the order of the model
%         RedOrder = order(Kmu1);
%        KmuR = balred(Kmu1,[RedOrder-1:-1:RedOrder-5]);
%%        
        KmuR = balred(Kmu1,13);%RedOrder
        figure;
        bodeplot(Kmu1,KmuR,P)
        legend('Full-order model','Reduced-order model');
        sysS = 1/(1+sys*KmuR);                
        sysSsample = usample(clp,1000);
        for k=1:1000
            if max(real(pole(sysSsample(:,:,k))))>=0
                display('unstable CL')
            end
        end        
%%         Plot closed-loop frequency response
        figure
        sysCLFreq =  sys*KmuR/(1+sys*KmuR); 
        systmp = usample(sysCLFreq,10);
        for kk = 1:10
            [magM,phaseM] = bode(systmp(:,:,kk,1),freqM*2*pi);
            magM = 20*log10(magM(:));
            phaseM = phaseM(:);
            for i= 1:size(phaseM,1)
                if phaseM(i,1)>400
                    phaseM(i,1) = phaseM(i,1)-720;
                elseif phaseM(i,1)>50
                    phaseM(i,1) = phaseM(i,1)-360;
                end
            end
%             subplot(2,1,1);
            semilogx(freqM,magM,'Linewidth',1);hold on;grid on;
            ylabel('Magnitude (dB)')
%             subplot(2,1,2);semilogx(freqM,phaseM,'Linewidth',1);hold on;grid on;
            xlabel('Frequency (Hz)')
%             ylabel('Phase (deg)')
        end
%%

        bode(Weinv,inv(Wu),sysS,KmuR*sysS)%Kmu,1-sysS,,sys*sysK/(1+sys*Kmu),P
        legend('Weighting function Inv','Weighting function u Inv','Senstivity','Control Input')%,,'Complementary sensitivity','Controller','Closed-loop freq resp'
        grid on;  
        
        %% Convert to discrete-time model for implementation on dSPACE
        KmuRD = c2d(KmuR,Ts,'prewarp',100*2*pi);
        
        BSfilterD = c2d(BSfilter,Ts,'prewarp',256*2*pi);
        [BSfilterD_num,BSfilterD_den] = tfdata(BSfilterD,'v');
%          BSfilterD = c2d(BSfilter,Ts,'matched');
%         BSfilterMod1D = c2d(BSfilterMod1,Ts,'prewarp',65*2*pi);      
        Kmu = KmuR*BSfilter; %
        % Determine whether to add the band stop filter 
        if BSFILTER == 1
            Kmu = KmuR * BSfilter;
            KmuD = KmuRD * BSfilterD;
        elseif BSFILTER == 0
            Kmu = KmuR;
            KmuD = KmuRD;
        end
        KmuD = balreal(KmuD);
        figure;
        bodeplot(Kmu,KmuD,P);
        legend('Continous SS','Discrete TF');
        [numKmuD,denKmuD] = tfdata(KmuD,'v');    
        [numKmu,denKmu] = tfdata(Kmu,'v'); 
        sysK = KmuD;
        end
end

%% Bode & Nyquist plot with model
figure;
[magK,phaseK] = bode(sysK,freqM*2*pi);
phaseK = phaseK(:);
for i= 1:size(phaseK,1)
    if phaseK(i,1)>20
        phaseK(i,1) = phaseK(i,1)-360;
    elseif phaseK(i,1)<-360
        phaseK(i,1) = phaseK(i,1)+360;
    end
end
magK = magK(:);
magTot = magM+20*log10(magK);
phaseTot = phaseM+phaseK(:);
for i= 1:size(phaseTot,1)
    if phaseTot(i,1)<-360
        phaseTot(i,1) = phaseTot(i,1)+360;
    elseif phaseTot(i,1)>100
        phaseTot(i,1) = phaseTot(i,1)-360;
    end
end
subplot(2,2,1);
semilogx(freqM,magTot,'r');
ylabel('Magnitued/dB'), hold on; grid on;
subplot(2,2,3);
semilogx(freqM,phaseTot,'r');
hold on;
semilogx(freqM,-180*ones(size(freqM)),'k');
semilogx(freqM,180*ones(size(freqM)),'k');
ylabel('Phase/deg')
xlabel('Frequency(Hz)'), grid on;

% nyquist plot
subplot(2,2,[2 4] )
Amp = 10.^(magTot/20);
j = sqrt(-1);
pha = phaseTot*pi/180;
freqresp = Amp.*exp(j*pha);
frdata = frd(freqresp,freqM*2*pi,0);
nyquistplot(frdata,Pn)
hold on;
scatter(-1,0,'ro')
axis auto
title('Nyqusit plot with Model');

return;
%% Bode & Nyquist plot with data
EXTRA_NQST =  0;
if EXTRA_NQST == 0
    if FREQ_NQST == 0
    freq = freq1;
    mag = magM;
    phase = phaseM;
    elseif FREQ_NQST == 1
        freq = freq1;
        mag = mag1;
        phase = phase1;
    elseif FREQ_NQST == 2
        freq = freq2;
        mag = mag2;
        phase = phase2;
    elseif FREQ_NQST == 3
        freq = freq3;
        mag = mag3;
        phase = phase3;
    elseif FREQ_NQST == 4
        freq = freq4;
        mag = mag4;
        phase = phase4;
    elseif FREQ_NQST == 5
        freq = freq5;
        mag = mag5;
        phase = phase5;
    end
elseif EXTRA_NQST == 2   
    load T2_100.X %T2_100
    freq3 = T2_100;
    load T2_100P.txt
    phase3 = T2_100P;
    for i= 1:size(phase3,1)
        if phase3(i,1)>20
            phase3(i,1) = phase3(i,1)-360;
        end
    end
    load T2_100M.txt
    mag3 = T2_100M-20*log10(100);%the output was multiplied by 100
    mag3 = mag3+20*log10(200*0.88);%the input was divided by 200 to get desired torque
    freq = freq3;
    mag = mag3;
    phase = phase3;
elseif EXTRA_NQST == 3
    load T3_1002.X 
    freq3 = T3_1002;
    load T3_1002P.txt
    phase3 = T3_1002P;
    for i= 1:size(phase3,1)
        if phase3(i,1)>20
            phase3(i,1) = phase3(i,1)-360;
        end
    end
    load T3_1002M.txt
    mag3 = T3_1002M-20*log10(100);%the output was multiplied by 100
    mag3 = mag3+20*log10(200);%the in
    freq = freq3;
    mag = mag3;
    phase = phase3;
elseif EXTRA_NQST == 4
    load T4_1002.X %T4_1002
    freq3 = T4_1002;
    load T4_1002P.txt
    phase3 = T4_1002P;
    for i= 1:size(phase3,1)
        if phase3(i,1)>20
            phase3(i,1) = phase3(i,1)-360;
        end
    end
    load T4_1002M.txt
    mag3 = T4_1002M-20*log10(100);%the output was multiplied by 100
    mag3 = mag3+20*log10(200);%the in
    freq = freq3;
    mag = mag3;
    phase = phase3;
elseif EXTRA_NQST == 5
    load T5_100.X %T5_100
    freq3 = T5_100;
    load T5_100P.txt
    phase3 = T5_100P;
    for i= 1:size(phase3,1)
        if phase3(i,1)>20
            phase3(i,1) = phase3(i,1)-360;
        end
    end
    load T5_100M.txt
    mag3 = T5_100M-20*log10(100);%the output was multiplied by 100
    mag3 = mag3+20*log10(200);%the in
    freq = freq3;
    mag = mag3;
    phase = phase3;
end

figure;
[magK,phaseK] = bode(sysK,freq*2*pi);
phaseK = phaseK(:);
for i= 1:size(phaseK,1)
    if phaseK(i,1)>20
        phaseK(i,1) = phaseK(i,1)-360;
    elseif phaseK(i,1)<-360
        phaseK(i,1) = phaseK(i,1)+360;
    end
end
magK = magK(:);
magTot = mag+20*log10(magK);
phaseTot = phase+phaseK(:);
for i= 1:size(phaseTot,1)
    if phaseTot(i,1)<-360
        phaseTot(i,1) = phaseTot(i,1)+360;
    elseif phaseTot(i,1)>100
        phaseTot(i,1) = phaseTot(i,1)-360;
    end
end

subplot(2,2,1);
semilogx(freq,magTot,'r');
ylabel('Magnitued/dB'), hold on; grid on;
subplot(2,2,3);
semilogx(freq,phaseTot,'r');
hold on;
semilogx(freq,-180*ones(size(freq)),'k');
semilogx(freq,180*ones(size(freq)),'k');
ylabel('Phase/deg')
xlabel('Frequency(Hz)'), grid on;

%%nyquist plot with data 
subplot(2,2,[2 4] )
Amp = 10.^(magTot/20);
j = sqrt(-1);
pha = phaseTot*pi/180;
freqresp = Amp.*exp(j*pha);
frdata = frd(freqresp,freq*2*pi,0);
nyquistplot(frdata,Pn)
hold on;
scatter(-1,0,'ro')
axis auto
title('Nyqusit plot with data');


% %band stop filter
% [magK,phaseK] = bode(BSfilterD,freq*2*pi);
% phaseK = phaseK(:);
% magK = magK(:);
% for i= 1:size(phaseK,1)
%     if phaseK(i,1)>20
%         phaseK(i,1) = phaseK(i,1)-360;
%     elseif phaseK(i,1)<-360
%         phaseK(i,1) = phaseK(i,1)+360;
%     end
% end
% subplot(2,1,1);
% semilogx(freq,20*log10(magK),'b');
% ylabel('Magnitued/dB'), grid on;hold off
% subplot(2,1,2);
% semilogx(freq,phaseK,'b');
% hold on;
% semilogx(freq,-180*ones(size(freq)),'k');
% semilogx(freq,180*ones(size(freq)),'k');
% ylabel('Phase/deg')
% xlabel('Frequency(Hz)'), grid on; hold off



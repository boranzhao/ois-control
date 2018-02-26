%*************************************************************************
%  
%  Classical controller design for miniaturized optical image stabilizers 
%             ----------------------------------------
%
% By Pan Zhao, Control Engineering Lab, UBC, BC, Canada
% Under supervision of Dr. Ryozo Nagamune.
% Revision: Feb 25, 2018
% Creation: Jan 21, 2016.

% Function: Main function for designing a classical controller consisting of a lead-lag compensator and three notch filters.
% Reference: Zhao, R. Nagamune and M. Chiao, \Multiple parameter-dependent robust control of miniaturized optical
% image stabilizers," accepted by Control Engineering Practice, 2018.

clc; clear; close all;

%% load frequency resposne data
% bode plot setting
P = bodeoptions; % Set phase visiblity to off and frequency units to Hz in options
P.FreqUnits = 'Hz'; % Create pl

for foldthecode=1
% LP 1
load LP1_100N.X %  LP1_100N
freq1 = LP1_100N;
load LP1_100NP.txt
phase1 = LP1_100NP;
for i= 1:size(phase1,1)
    if phase1(i,1)>20
        phase1(i,1) = phase1(i,1)-360;
    end
end
load LP1_100NM.txt
mag1 = LP1_100NM-20*log10(100);%the output was multiplied by 100
mag1 = mag1 + 20*log10(200*1.4287);%the input was divided by 200 to get desired torque
subplot(2,1,1);
semilogx(freq1,mag1,'Linewidth',2);hold on;
grid on;
subplot(2,1,2);
semilogx(freq1,phase1,'Linewidth',2);hold on;
grid on;

% LP 2
load LP2_100N.X
freq2 = LP2_100N;
load LP2_100NP.txt
phase2 = LP2_100NP;
for i= 1:size(phase2,1)
    if phase2(i,1)>20
        phase2(i,1) = phase2(i,1)-360;
    end
end
load LP2_100NM.txt
mag2 = LP2_100NM-20*log10(100);%the output was multiplied by 100
mag2 = mag2+ 20*log10(200);%the input was divided by 200 to get desired torque
% mag = mag+20*log10(400);%Use a gain of 700
subplot(2,1,1);
semilogx(freq2,mag2,':g','Linewidth',2);hold on;
grid on;
subplot(2,1,2);
semilogx(freq2,phase2,':g','Linewidth',2);hold on;
grid on;

% LP3
load LP3_100.X %LP3_100
freq3 = LP3_100;
load LP3_100P.txt
phase3 = LP3_100P;
for i= 1:size(phase3,1)
    if phase3(i,1)>20
        phase3(i,1) = phase3(i,1)-360;
    end
end
load LP3_100M.txt
mag3 = LP3_100M-20*log10(100);%the output was multiplied by 100
mag3 = mag3+20*log10(200*0.88);%the input was divided by 200 to get desired torque
subplot(2,1,1);
semilogx(freq3,mag3,'r','Linewidth',2,'MarkerSize',2);
ylabel('Magnitued (dB)');
grid on;
subplot(2,1,2);
semilogx(freq3,phase3,'r','Linewidth',2,'MarkerSize',2);
ylabel('Phase (deg)')
grid on;

% LP 4
load LP4_100N.X
freq4 = LP4_100N;
load LP4_100NP.txt
phase4 = LP4_100NP;
for i= 1:size(phase4,1)
    if phase4(i,1)>10
        phase4(i,1) = phase4(i,1)-360;
    end
end
load LP4_100NM.txt
mag4 = LP4_100NM-20*log10(100);%the output was multiplied by 100
mag4 = mag4+20*log10(200);%the input was divided by 200 to get desired torque
subplot(2,1,1);
semilogx(freq4,mag4,':c','Linewidth',2,'MarkerSize',2);
ylabel('Magnitued (dB)');
grid on;
subplot(2,1,2);
semilogx(freq4,phase4,':c','Linewidth',2,'MarkerSize',2);
ylabel('Phase (deg)')
grid on;
subplot(2,1,1);
xlabel('Frequency (Hz)','Fontsize',13);
legend('Data1','Data2','Data3','Data4')
subplot(2,1,2);
xlabel('Frequency (Hz)','Fontsize',13);
legend('Data1','Data2','Data3','Data4')

% LP5
load LP5_1002N.X
freq5 = LP5_1002N;
load LP5_1002NP.txt
phase5 = LP5_1002NP;
for i= 1:size(phase5,1)
    if phase5(i,1)>10
        phase5(i,1) = phase5(i,1)-360;
    end
end
load LP5_1002NM.txt
mag5 = LP5_1002NM-20*log10(100);%the output was multiplied by 100
mag5 = mag5+20*log10(200);%the input was divided by 200 to get desired torque
subplot(2,1,1);
semilogx(freq5,mag5,'c','Linewidth',2,'MarkerSize',2);
ylabel('Magnitued (dB)');
grid on;
subplot(2,1,2);
semilogx(freq5,phase5,'c','Linewidth',2,'MarkerSize',2);
ylabel('Phase (deg)')
grid on;
subplot(2,1,1);
xlabel('Frequency (Hz)','Fontsize',13);
legend('Data1','Data2','Data3','Data4','Data5')
subplot(2,1,2);
xlabel('Frequency (Hz)','Fontsize',13);
legend('Data1','Data2','Data3','Data4','Data5')
end

%% Choose which data for validation
freq = freq1*2*pi;
mag = mag1;
phase = phase1;

%% Controller Design 
% add the notch filters
j= sqrt(-1);
% good control parameter for nomial plant
w = [62.84*2*pi 97*2*pi 100*2*pi 105*2*pi 250*2*pi 271*2*pi]; %1592 664.4 1658.1 1660  
alpha = [20 5 5 5 20 20];%100

Knf12= tf(1);
for i=1:length(w)-2    
    beta(i) = sqrt(alpha(i)^2+w(i)^2);
    modeltemp = tf(poly([-alpha(i)+w(i)*j -alpha(i)-w(i)*j]), poly([-beta(i) -beta(i)]));
    Knf12 =  Knf12*modeltemp;
end
Knf3 = tf(1);
for i=length(w)-1:length(w)   
    beta(i) = sqrt(alpha(i)^2+w(i)^2);
    modeltemp = tf(poly([-alpha(i)+w(i)*j -alpha(i)-w(i)*j]), poly([-beta(i) -beta(i)]));
    Knf3 =  Knf3*modeltemp;
end
Knf = Knf12*Knf3;
% Knf = Knf*tf(1000^3, poly([-1000 -1000 -1000]));
[magN,phaseN] = bode(Knf,freq);
phaseN = phaseN(:);
magN = magN(:);
% 
% [magN1,phaseN1] = bode(Model,freq);
% phaseN1 = phaseN1(:);
% magN1 = magN1(:);
% magN = magN.*magN1;
% phaseN = phaseN+phaseN1;
for i= 1:size(phaseN,1)
    if phaseN(i,1)>20
        phaseN(i,1) = phaseN(i,1)-360;
    end
end
subplot(2,1,1);
semilogx(freq/2/pi,20*log10(magN)+mag,'g');%+mag
ylabel('Magnitued/dB'), grid on;hold on;
subplot(2,1,2);
semilogx(freq/2/pi,phaseN(:)+phase,'g');%+phase
ylabel('Phase/deg')
xlabel('Frequency(Hz)'), grid on; hold on;

% add the 1st lag compensator
% sysC = tf([1/21.59 1],[1/6.14 1])*2; %lag compensator
% sysC = tf([1/51.59 1],[1/6.14 1])*2; %lag compensator
Klag = tf([1/150.59 1],[1/25 1])*1.7; %lag compensator
% Klag= tf(1);
[magC,phaseC] = bode(Klag,freq);
phaseC = phaseC(:);
magC = magC(:);
for i= 1:size(phaseC,1)
    if phaseC(i,1)>20
        phaseC(i,1) = phaseC(i,1)-360;
    end
end
magTot = 20*log10(magC)+20*log10(magN)+mag;
phaseTot = phaseC(:)+phaseN(:)+phase;
subplot(2,1,1);
semilogx(freq/2/pi,magTot,'m');
ylabel('Magnitued/dB'), grid on;
subplot(2,1,2);
semilogx(freq/2/pi,phaseTot,'m');
ylabel('Phase/deg')
xlabel('Frequency(Hz)'), grid on; 
semilogx(freq/2/pi,-180*ones(size(phase)),'k');
% 
% add the lead compensator
zl= 200.8;
pl = 1200;
% zl= 30.8;
% pl = 1200;
Klead = tf([1/zl 1],[1/pl 1]);
% Klead = tf(1);
[magL,phaseL] = bode(Klead,freq);
phaseL = phaseL(:);
magL = magL(:);
for i= 1:size(phaseL,1)
    if phaseL(i,1)>20
        phaseL(i,1) = phaseL(i,1)-360;
    elseif phaseL(i,1)<-360
        phaseL(i,1) = phaseL(i,1)+360;
    end
end
magTot = magTot+20*log10(magL);
phaseTot = phaseTot+phaseL(:);
for i= 1:size(phaseTot,1)
    if phaseTot(i,1)<-360
        phaseTot(i,1) = phaseTot(i,1)+360;
    end
end
subplot(2,1,1);
semilogx(freq/2/pi,magTot,'r');
ylabel('Magnitued/dB'), grid on;hold off
legend('origin.','origin.+notch.','origin.++notch+lag.','origin.+notch+lag.+lead')
subplot(2,1,2);
semilogx(freq/2/pi,phaseTot,'r');
hold on;
semilogx(freq/2/pi,-180*ones(size(phase)),'k');
ylabel('Phase/deg')
xlabel('Frequency(Hz)'), grid on; hold off
legend('origin.','origin.+notch.','origin.++notch+lag.','origin.+notch+lag.+lead')
legend;

%% Nyquist plot to verify the stability
Amp = 10.^(magTot/20);
j = sqrt(-1);
pha = phaseTot*pi/180;
freqresp = Amp.*exp(j*pha);
frdata = idfrd(freqresp,freq,0);
figure;
nyquist(frdata)
hold on;
scatter(-1,0,'ro')
axis auto

%% Convert to discrete-time model for implementation on dSPACE
Kcla= Knf12*Knf3*Klag*Klead;
Kcla = ss(Kcla);
Kcla = balreal(Kcla);


Ts = 5e-4; %simulation step
Knf12D = c2d(Knf12*Klag*Klead,Ts,'prewarp',101*2*pi);
Knf3D = c2d(Knf3,Ts,'prewarp',262*2*pi);
KclaD = Knf12D * Knf3D;
figure;
bodeplot(Kcla,KclaD,P)
legend('Continous','tustin')
[numKllRD,denKllRD] = tfdata(KclaD,'v');




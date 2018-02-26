%% Load and Plot frequency response data for all the plants

function ReadFreqResp()

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
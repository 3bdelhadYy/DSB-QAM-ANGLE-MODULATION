clear all; clc;

%Sampling
f_s=100000;
t_s=1/f_s;
t=0:t_s:0.002;

%carriar signal 
A_c=1;
f_c=10000;

%Modulation constant
K_p1=0.05;
K_p2=1;
K_p3=5;
K_p4=10;
K_f1=1000;
K_f2=2000;

%=====================(1)=====================

%Information signal m_1(t) --> Triangular signal
% --> use sawtooth() function to generate it
T_1=0.001;
f_1=1/T_1;
m_1=sawtooth(-2*pi*f_1*t + pi,1);

% %Information signal m_2(t) --> step signal 
m_2_Amp=[1,0.5,-0.5,-1,1];
m_2_steps=[0,0.0005,0.001,0.0015,0.002];
m_2=interp1(m_2_steps, m_2_Amp, t, 'previous');

%plots
figure;
subplot(2,1,1);
plot(t*1e3, m_1);
title('Information signal m_1(t)'); xlabel('Time (msec)'); ylabel('m_1(t)');
grid on;

subplot(2,1,2);
plot(t*1e3, m_2);
title('Information signal m_2(t)'); xlabel('Time (msec)'); ylabel('m_2(t)');
grid on;
%savefig('message_signals.fig');

%=====================(2)=====================

%Modulated signal s_1(t)
s_1=A_c* cos(2*pi*f_c*t+K_p1*m_1);

%plots
figure;
plot(t*1e3, s_1);
title('Modulated signal s_1(t)'); xlabel('Time (msec)'); ylabel('s_1(t)');
grid on
%savefig('modulated_signal.fig');
%=====================(3)=====================
s_1_2=A_c* cos(2*pi*f_c*t+K_p2*m_1);
s_1_3=A_c* cos(2*pi*f_c*t+K_p3*m_1);
s_1_4=A_c* cos(2*pi*f_c*t+K_p4*m_1);

%plots
figure;
subplot(3,1,1)
plot(t*1e3, s_1_2);
title('Modulated signal s_1(t) with K_p=1'); xlabel('Time (msec)'); ylabel('s_1(t)');
grid on

subplot(3,1,2)
plot(t*1e3, s_1_3);
title('Modulated signal s_1(t) with K_p=5'); xlabel('Time (msec)'); ylabel('s_1(t)');
grid on

subplot(3,1,3)
plot(t*1e3, s_1_4);
title('Modulated signal s_1(t) with K_p=10'); xlabel('Time (msec)'); ylabel('s_1(t)');
grid on
%savefig('modulated_signal_with_diff_Kp.fig');
%=====================(4)=====================

s_2=A_c*cos(2*pi*f_c*t+K_f1*cumsum(m_2)*t_s);
%s_2=A_c*cos(2*pi*f_c*t+K_f1*integral(m_2,inf,t));


figure;
plot(t*1e3, s_2);
title('Modulated signal s_2(t) with K_f=1000'); xlabel('Time (msec)'); ylabel('s_2(t)');
grid on
%savefig('modulated_signal_with_Kf_1000.fig');
%=====================(5)=====================
s_3=A_c*cos(2*pi*f_c*t+K_f2*cumsum(m_1)*t_s);

figure;
plot(t*1e3, s_3);
title('Modulated signal s_3(t) with K_f=2000'); xlabel('Time (msec)'); ylabel('s_3(t)');
grid on
%savefig('modulated_signal_with_Kf_2000.fig');

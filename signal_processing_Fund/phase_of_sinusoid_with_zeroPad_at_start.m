clc;
clear all;
close all;

%
f1 = 100;
fs = 20000;
Ts = 1/fs;
T1 =  1/f1;
N_cycles = 100;
%t = 0:Ts:((N_cycles*T1)-213*Ts);
t1 = 0:Ts:((N_cycles*T1));


Tau = 0.90*T1;



t_pad_index = 0:Ts:(Tau);
t_pad_len = length(t_pad_index);



x = zeros(1,(t_pad_len +length(t1)) );

phase1 = (pi/4)*1;
%phase1 = (2*pi*f1*Tau);

x1 = cos(2*pi*f1*t1 + phase1 ); % It is wise to use cosine signal - Phase of DFT is calc wrt zero-phase cosine 

x(t_pad_len+1:end) = x(t_pad_len+1:end) + x1;

%x = x1;

time_delay_phase_shift = -(2*pi*f1*Tau);
signal_phase = phase1;
Expected_phase = time_delay_phase_shift + signal_phase;
fprintf("%f %f Expected_phase=%f\n",rad_2_deg(time_delay_phase_shift),rad_2_deg(signal_phase),rad_2_deg(Expected_phase));

time_delay_phase_shift_wr = mod(time_delay_phase_shift,(2*pi));
signal_phase_wr = mod(signal_phase,(2*pi));
Expected_phase_wr = time_delay_phase_shift_wr + signal_phase_wr;

fprintf("%f %f Expected_phase_wr=%f\n",rad_2_deg(time_delay_phase_shift_wr),rad_2_deg(signal_phase_wr),rad_2_deg(Expected_phase_wr));


t = 0:Ts:(length(x)-1)*Ts;
figure;
plot(t,fftshift(x));
%plot(t,x);
xlim([0,1.02]);
ylim([-2,2]);
%plot(t(1:550),x(1:550));
%xlim([0 550*Ts]);
%pause;
%display(mean(x)); --> notice here mean(x) does not correspond to dc comp wrt fft
                      % fft correlated samples over several cycles to evaluate freq.




FFT_SIZE = length(x) + mod(length(x),2);
FFT_SIZE_BY_2 = (FFT_SIZE/2);

fft_res = fs/FFT_SIZE;
k = 1:1:FFT_SIZE; % k starting from 1 is needed for referencing fft_op_array
f = fft_res*(k-1);

X = fft(x,FFT_SIZE);
%X = fft(fftshift(x),FFT_SIZE);
X_abs = abs(X);
X_angle = angle(X);


figure;
stem(f,X_abs);

figure;
fprintf("X_abs(1) =%d  %1.25f\n",X_abs(1),X_abs(1));
stem(f,X_abs);
xlim([-1,20]);
ylim([0 4.8]);


figure;
fprintf("X_abs(1) =%d  %1.25f\n",X_abs(1),X_abs(1));
stem(f,X_abs);
xlim([-1,f(8)]);
ylim([0 0.1]);


figure;
fprintf("X_abs(1) =%d  %1.25f\n",X_abs(1),X_abs(1));
stem(f,X_abs);
xlim([-1,f(3)]);
ylim([0 5e-13]);




figure;
stem(f,X_angle);


figure;
stem(f,X_angle);
xlim([-1,8]);
ylim([-pi pi]);

X_angle_deg = rad_2_deg(X_angle);
figure;
stem(f,X_angle_deg);
xlim([-1,80]);
ylim([-180 180]);

%angle(X)

MPH = 600;

[kpks1,klocs1] = findpeaks(X_abs(1:FFT_SIZE_BY_2) ,k(1:FFT_SIZE_BY_2) ,'MinPeakHeight',MPH);
display(kpks1);
display(klocs1);
display(klocs1*fft_res);

fprintf("phase of zero freq =%d  %2.20f\n",X_angle(1),X_angle(1));
display(X_angle(1));
display(rad_2_deg(X_angle(1)));

% here DC comp != 0 because of Rounding Error in FFT calculation.
% This Error is +ve decimal number so we get phase at f=0 to be 0.

display(X(1));
fprintf("fft of DC comp =%d  %1.25f\n",X(1),X(1));
fprintf("real of fft-DC comp =%d  %1.25f\n",real(X(1)),real(X(1)));
fprintf("imag of fft-DC comp =%d  %1.25f\n",imag(X(1)),imag(X(1)));
b_by_a = imag(X(1))/real(X(1))
atan( b_by_a )
atan2( imag(X(1)) , real(X(1)) )


display(X(klocs1));
fprintf("real of X(klocs1) =%d  %1.25f\n",real(X(klocs1)),real(X(klocs1)));
fprintf("imag of X(klocs1) =%d  %1.25f\n",imag(X(klocs1)),imag(X(klocs1)));
display(X_angle(klocs1));
theta_deg = rad_2_deg( X_angle(klocs1) );
display(theta_deg);
%Observe Here theta_deg= -85 unlike -90 in previous case .
%This is because of Rounding Error in fft X(100)=(3.7656 -48.3441i)
%is not accurate .

b_by_a = imag(X(klocs1))/real(X(klocs1))
phas = atan(imag(X(klocs1))/real(X(klocs1)))
rad_2_deg(phas)
phas = atan2(imag(X(klocs1)),real(X(klocs1)))
rad_2_deg(phas)

fprintf("@ f=klocs1 phase =%d  %2.20f deg=%2.20f\n",X_angle(klocs1),X_angle(klocs1),rad_2_deg(X_angle(klocs1)));


X_roundingNoiseFiltered = X;
indices_to_be_removed = 1:length(X);
indices_to_be_removed(klocs1) = [];
X_roundingNoiseFiltered(indices_to_be_removed) = 0;
X1_abs = abs(X_roundingNoiseFiltered);
X1_angle = angle(X_roundingNoiseFiltered);


figure,stem(f,X1_abs);

figure,stem(f,X1_abs);
ylim([0,1e-14]);

figure;
stem(f,rad_2_deg(X1_angle));
ylim([-180,180]);

%display(Expected_phase);
display(rad_2_deg(Expected_phase_wr));

return;

function x_deg = rad_2_deg(x_rad)
    x_deg = (x_rad)*(180/pi);
end

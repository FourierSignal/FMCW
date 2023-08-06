%%

% fs_adc  =  50MHZ;
% Tc = 100us => 5000 samples /chirp

% N = 4096  =>  fres = fs_adc/N = 12.4Khz
% RFFT's IF freq range ( 0 ----> 25Mhz)
%           dist range ( 0 ----> 5000 mtrs)


% But we fixed f_IF_detectable = 10Mhz  (as per our design actual fs = 25Mhz)
% That means maximum detectable range is 2000 mtrs
% so do not include sinusoids with f_Ifs > 10Mhz / range > 2000 mtrs in the simulation. 



%%


clc;
clear all;
close all;

%%
GHz = 1e9;
MHz = 1e6;
us = 1e-6;
KHz = 1e3;
C = 3e8;


%%
%This program generates IF Signal(of duration one Frame) Sampled(ADC output)
%at 50MHz rate.
% It simulates Reflections from moving Target @Range R , with velocity v.
%%

f_centre = 16.5*GHz;
f_sweep =  75*MHz;
f_low  = f_centre - (f_sweep/2);
f_high = f_centre + (f_sweep/2);

lambda = C/f_centre;
%N_beams = 23;
N_beams = 1;
T_chirp = 100*us;
Sweep_rate = f_sweep/T_chirp;
fprintf("Sweep_rate=%f Hz/sec %f MHz/usec\n",Sweep_rate, Sweep_rate/MHz  );
N_chirps = 512;
%N_chirps = 8;


%
RFFT = 4096; % to achieve FFT Range Resol = 2.4414 mtrs
RFFT_BY_2 = RFFT/2;
DFFT = N_chirps;
%

fs = 25*MHz;
f_IF_instr = fs/2.5;

R_instr =  f_IF_instr * (C/(2*Sweep_rate));
fprintf("f_IF_instr=%d   R_instr=%d\n",f_IF_instr,R_instr);

R_det_max =  R_instr - (0.1 * R_instr);
f_IF_det = ((2*Sweep_rate)/C) * R_det_max;

fprintf("f_IF_det=%d   R_det_max=%d\n",f_IF_det,R_det_max); 
fprintf("Enter Range within R_det_max=%d\n",R_det_max);


%sim_freq_array = input("enter array of Freqs to be simulated >");
sim_range_array = input("enter array of Ranges of targets to be simulated >");
sim_vel_array = input("enter array of velocities of these targets >");
display(sim_range_array);


fs_adc = 50*MHz;
Ts = 1/fs_adc;
T_frame = (N_chirps * T_chirp);
Frame_size = (T_frame/Ts);
N_max  = uint32(T_chirp/Ts); % max number of samples available per chirp

A = 10*(10^-3);


IF_Sig_SIZE = uint32(N_beams*N_chirps*N_max);
Simulated_IF_Sig_double_Nmax = zeros(1,IF_Sig_SIZE); % for 5000 samples

Simulated_IF_Sig_double_Matrix = zeros(DFFT,RFFT);



index_t = 0 : Ts : T_chirp-Ts;

for beam = 1 : N_beams   
    for chirp = 1 : N_chirps 
        Simulated_chirp_IFSig_Nmax_double = zeros(1,uint32(N_max));
        for i = 1:length(sim_range_array)
            display(chirp);
            %display(i);
            
            range = sim_range_array(i);
            Tau = (2*range)/C;
            f_IF = 2*Sweep_rate*(range/C);
            Phase_IF = 4*pi*(range/lambda); 
            velocity = sim_vel_array(i);
            delta_Tau = (2*(velocity*T_chirp))/C;
            delay_Tau = Tau + ((chirp-1)*delta_Tau);

            fprintf("range=%d Tau=%d velocity=%d delay_Tau=%d\n",range,Tau,velocity,delay_Tau);

            delta_IF = ((2*Sweep_rate)/C)*(velocity*T_chirp);
            delta_phase = ((4*pi)/lambda)*(velocity*T_chirp);
           
            %fprintf("range=%f Tau=%f\n",range,Tau);
            %fprintf("f_IF=%f delta_IF=%f\n",f_IF,delta_IF);
            fprintf("Phase_IF=%f delta_phase=%f\n",Phase_IF,delta_phase);
            fprintf("Phase_IF=%f delta_phase=%f\n",rad_2_deg(Phase_IF),rad_2_deg(delta_phase));
        
            f =  f_IF + (chirp-1)* delta_IF;

            %phase_angle =  ceil(Phase_IF + ((chirp-1)*delta_phase));
            %This is a mistake ...that change in phase caused by ceil is
            %creating several peaks in dopplerFFT.
            %observed that This causes steps in Phase response which will
            %be removed when ceil is removed.
            % study about this Phase change effects on dopplerFFT output.
            phase_angle =  Phase_IF + ((chirp-1)*delta_phase);
            wrapped_phase_angle = mod(ceil(phase_angle) ,(2*pi));
            %fprintf("f=%f phase_angle=%f wrapped_phase_angle=%f\n",f,phase_angle,wrapped_phase_angle);
fprintf("f=%f phase_angle=%f wrapped_phase_angle=%f\n",f,rad_2_deg(phase_angle),rad_2_deg(wrapped_phase_angle));
            
            
            %index_t = Tau : Ts : T_chirp-Ts;
            Simulated_chirp_IFComp_Sig_Nmax_double = zeros(1,N_max);

            ifComp_start = floor(delay_Tau/Ts);
            ifComp_Timeslot = T_chirp - delay_Tau;
           
            index_t_ifComp = 0 : Ts : ifComp_Timeslot;
            Simulated_chirp_IFComp_Sig_Nmax_double(ifComp_start+1:N_max) = A * cos((2 * pi * f * index_t_ifComp) + phase_angle);
            
            
            
            %{
            figure;
            ifComp_start = floor(delay_Tau/Ts)-10;
            ifComp_end = ifComp_start + 40;
            plot(index_t(ifComp_start:ifComp_end),Simulated_chirp_IFComp_Sig_Nmax_double(ifComp_start:ifComp_end));
            pause;
            %}

            %N_IfComp = length(index_t);
            %start_ifComp = (N_max - N_IfComp)+1;
            %fprintf("N_IfComp=%d start_ifComp=%d N_max=%d\n",N_IfComp,start_ifComp,N_max);
            %cosine_index_t = 0 : Ts : ((N_IfComp-1)*Ts);
            %Simulated_chirp_IFComp_Sig_Nmax_double(start_ifComp:N_max) = A * cos((2 * pi * f * cosine_index_t));
            Simulated_chirp_IFSig_Nmax_double = Simulated_chirp_IFSig_Nmax_double + Simulated_chirp_IFComp_Sig_Nmax_double;
            %
        end
        
        start_n =  uint32(((chirp-1)*N_max)+1);
        last_n =   uint32(chirp*N_max);
        Simulated_IF_Sig_double_Nmax(start_n : last_n) = Simulated_chirp_IFSig_Nmax_double;

        Simulated_IF_Sig_double_Matrix(chirp, : ) = Simulated_chirp_IFSig_Nmax_double(1:RFFT);


     end
    
end

Simulated_IF_Sig_double = reshape(Simulated_IF_Sig_double_Matrix',[1,RFFT*DFFT]);

figure;
%plot(Simulated_IF_Sig_double_Nmax(1:400));
%plot(index_t,Simulated_IF_Sig_double_Nmax(1:N_max));
stem(Simulated_IF_Sig_double_Nmax(1:3*N_max));

str1 = '';
for i = 1: length(sim_range_array)
   str2 = sprintf("%d_v_%d_",sim_range_array(i),sim_vel_array(i));
    %display(str2); 
    str1 = strcat(str1, str2);
end
file_name = strcat('Simulated_Sig_',str1,'double.bin')
%return;

rel_dir = "H:\DSP_related\MatLab_Ra_design_structured";
ip_file_dir = fullfile(rel_dir,"\floatingPt_Sim_adc_data");

ip_file3 =  fullfile(ip_file_dir,file_name);
disp(ip_file3);


[fileID3,errmsg] = fopen(ip_file3, 'w');
if fileID3 < 0
    disp(errmsg);
    return;
end

bytes_written = fwrite(fileID3,Simulated_IF_Sig_double,"double");
fclose(fileID3);
fprintf("bytes_written=%d\n",bytes_written);

return;


function x_deg = rad_2_deg(x_rad)
    x_deg = (x_rad)*(180/pi);
end

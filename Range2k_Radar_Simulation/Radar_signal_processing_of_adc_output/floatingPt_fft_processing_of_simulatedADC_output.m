% 
% floating point processing of simulated adc data 
% removed windowing and lpf 
% RFFT = 4096
% fs_adc = 50Mhz

%%
clear;
close all;
clc;

%%
MHz = 1e6;
us = 1e-6;
KHz = 1e3;
C = 3e8;
%% Chirp Signals
fc = 16.5e9;
lambda = C/fc;
f_sweep = 75*MHz;
sweep_time = 100*us;
T_chirp = 100*us; % some times this can be slightly more than sweep_time
N_chirps = 512;
Sweep_rate = f_sweep/sweep_time;
fprintf("Sweep_rate=%f Hz/sec %f MHz/usec\n",Sweep_rate, Sweep_rate/MHz  );

fs = 25*MHz;
ts = 1/fs;
%%
RFFT = 4096;
RFFT_BY_2 = RFFT/2;
DFFT = N_chirps;
%%
fs_adc = 50*MHz;
Ts_adc = 1/fs_adc;

R_det_max = 1.8e3;
f_IF_det = ((2*Sweep_rate)/C) * R_det_max;
fprintf("f_IF_det=%d   R_det_max=%d\n",f_IF_det,R_det_max); 

R_instr =  R_det_max + (0.1 * R_det_max);
f_IF_instr =  ((2*Sweep_rate)/C) * R_instr;
fprintf("f_IF_instr=%d   R_instr=%d\n",f_IF_instr,R_instr);

R_instr = 2e3;

%%

%[bm,ak] = define_window_based_fir_filter(fs);

rel_dir = "H:\DSP_related\MatLab_Ra_design_structured";
ip_file_dir = fullfile(rel_dir,"\floatingPt_Sim_adc_data");
%ip_filename =  fullfile(ip_file_dir,'Simulated_Sig1_700_double.bin');

ip_filename =  fullfile(ip_file_dir,'Simulated_Sig_700_v_8_double.bin');

%ip_filename =  fullfile(ip_file_dir,'Simulated_Sig_700_v_8_800_v_35_double.bin');


disp(ip_filename);
[fid,errmsg] = fopen(ip_filename, 'rb');
if fid < 0
    disp(errmsg);
    return;
end

adc_out = fread(fid, 'double');
adc_out = adc_out';
fclose(fid);


%adc_out = downsample(x_data,fix((2*fs)/fs)); % needed for non simulated Sigals. 
%adc_out = test_data_buf;


%{
scf = input('Provide the frame number for the ADC data to be analyzed   ');
L = RFFT_BY_2;
y = adc_out(((scf-1)*L+1):(scf+2)*L);
figure(1);
n = 1:1:L;   
ts = 1/fs;
t = n*ts;
stem(y(1:400),'r','LineWidth',1);
%stem(y,'r','LineWidth',1);
str = ['Samples of ADC data frame of data frame number ',num2str(scf)];
title(str, 'FontWeight','Bold','FontSize',16);
xlabel('Time Domain in sec', 'FontWeight','Bold','FontSize',12)
return;
%}

N_max  = uint32(T_chirp/Ts_adc);

frame_size = RFFT;
beam_size = DFFT*frame_size;
NUM_FRAMES= fix(length(adc_out)/frame_size);
NUM_BEAMS =fix( length(adc_out)/(beam_size));
fprintf("NUM_BEAMS= %d NUM_FRAMES= %d\n",NUM_BEAMS,NUM_FRAMES);

k_range= 1:1:RFFT;
fres_rfft = fs_adc/RFFT; 
freq_vect = k_range*fres_rfft;
r_res_rfft = fres_rfft * ((C*T_chirp)/(2*f_sweep));
range_vect = k_range * r_res_rfft;

win_rfft = window(@hamming, RFFT);
for beam_num = 1:1:NUM_BEAMS
    x_beam_data = adc_out(((beam_num-1)*beam_size+1):(beam_num*beam_size));
    for chirp_index = 1:1:DFFT
        display(chirp_index);
        
        %x_seg = x_beam_data((chirp_index-1)*N_max+1:chirp_index*N_max);
        x_seg = x_beam_data((chirp_index-1)*frame_size+1:chirp_index*frame_size);
        
        x_chirp = x_seg( 1 : uint32(frame_size) );    
       
        % needed  when Signal-duration is non multiple of it's period.
        % x_filt = filter(bm, ak,x_chirp);
        % y = x_filt .* win_rfft;
        
        y = x_chirp;
        rangefft_op_of_chirp = fft(y, RFFT);
        
        rangefft_op_of_chirp_abs = abs(rangefft_op_of_chirp(1:RFFT_BY_2));
        MPH = 12;
        [kpks1,klocs1] = findpeaks(rangefft_op_of_chirp_abs,range_vect(1:RFFT_BY_2),'MinPeakHeight',MPH);
        display(klocs1);
        
        %{
        figure;
        findpeaks(rangefft_op_of_chirp_abs,range_vect(1:RFFT_BY_2),'MinPeakHeight',MPH)
        pause;
        %}

        Rangefft_op_matrix(chirp_index, :) = rangefft_op_of_chirp;
        Rangefft_matrix(chirp_index, :) = rangefft_op_of_chirp(1:RFFT_BY_2);

    end
end




Rfft_op_avg_vector =  mean(abs(Rangefft_op_matrix));
figure;
MPH = 2;
findpeaks(Rfft_op_avg_vector,'MinPeakHeight',MPH);
xlim([0 RFFT]);
ylim([0 20]);


range_det_str_vect = mean(abs(Rangefft_matrix), 1); % mean along-dim1 / mean of each column
k_range= 1:1:RFFT;
[kpks,klocs] = findpeaks(range_det_str_vect,k_range(1:RFFT_BY_2),'MinPeakHeight',MPH);
%display(locs);



% use Instrumented Range for Plotting 

fres_rfft = fs_adc/RFFT; 
freq_vect = k_range*fres_rfft;
r_res_rfft = fres_rfft * ((C*T_chirp)/(2*f_sweep));
range_vect = k_range * r_res_rfft;

k_instr = ceil(R_instr / r_res_rfft);
%k_instr = k_range(k_max);
display(k_instr);

%{
figure;
plot(k_range(2:k_instr),range_det_str_vect(2:k_instr));
title('Range Detectiopn Strength Vs k ', 'FontWeight','Bold','FontSize',16);
xlabel('k', 'FontWeight','Bold','FontSize',12)
ylabel('Magnitude Spectrum values', 'FontWeight','Bold','FontSize',12);
xlim([0 k_instr]);
ylim([0 20]);


figure;
plot(freq_vect(2:k_instr),range_det_str_vect(2:k_instr));
xlim([0 freq_vect(k_instr)]);
ylim([0 20]);
%}

figure;
plot(range_vect(2:k_instr),range_det_str_vect(2:k_instr));
xlim([0 range_vect(k_instr)]);
ylim([0 20]);




Rangefft_matrix_transpose = Rangefft_matrix';

rfft_of_rcell = Rangefft_matrix_transpose(klocs(1),:);

chirp_index = 1 : 1 : DFFT;
figure;
stem (chirp_index,abs(rfft_of_rcell));
%xlim([1,10]);

figure;
stem (chirp_index,angle(rfft_of_rcell));
%xlim([1,10]);

figure;
plot (chirp_index,angle(rfft_of_rcell));
%xlim([1,10]);

figure;
plot (chirp_index,angle(rfft_of_rcell));
xlim([1,50]);

phase_signal = angle(rfft_of_rcell);



fs_doppler = 1/sweep_time; % every beat is sampled at chirp repetetion freq.
doppler_fft_res = fs_doppler/DFFT; % doppler_fft_res = 1/(Tc*Nc)
doppler_freq =  0 : doppler_fft_res : (DFFT -1)*(doppler_fft_res); 

doppler_index = 1 : 1 : DFFT;
doppler_vel_res = (lambda/2)*doppler_fft_res;
doppler_vel = doppler_index * doppler_vel_res; 

doppler_index_symm = doppler_index - (DFFT/2) ;
doppler_vel_symm = doppler_index_symm * doppler_vel_res; 

k_dopp = doppler_index_symm;


h_d = chebwin(DFFT);
dopplerfft_matrix = zeros(RFFT_BY_2,DFFT);
for rcell = 1: RFFT_BY_2
    kth_range_rfft_samples_over_frame = Rangefft_matrix_transpose(rcell,:) .* h_d';
    rcell_dfft = fft(kth_range_rfft_samples_over_frame, DFFT);
    %rcell_dfft = fft(angle(kth_range_rfft_samples_over_frame), DFFT);
    dopplerfft_matrix(rcell, :) = fftshift(rcell_dfft);

    if ismember(rcell,klocs) %|| rcell == 10 || rcell < 6
        %{
        figure;
        plot(abs(Rangefft_matrix_transpose(rcell,:)));
        title(['rcell = ',num2str(rcell)],'Fontsize',16);

        figure;
        plot(angle(Rangefft_matrix_transpose(rcell,:)));
        title(['rcell = ',num2str(rcell)],'Fontsize',16);
        pause;

        %}

        %
        rcell_dfft_abs = abs(dopplerfft_matrix(rcell,:));

        max_val = max(rcell_dfft_abs);
        index = find(rcell_dfft_abs == max_val);
        if  max_val > 1000 
            fprintf("rcell=%d  max=%d index=%d\n",rcell,max_val,index);
        end

        figure;
        plot(doppler_vel_symm,abs(dopplerfft_matrix(rcell,:)));
        ylim([0,max_val+10]);
        %ylim([0,0.15]);  % Why there are spikes in range cells where there is no target ??
        title(['rcell = ',num2str(rcell)],'Fontsize',16);

        figure;
        plot(doppler_vel_symm,abs(dopplerfft_matrix(rcell,:)));
        %ylim([0,max_val+10]);
        ylim([0,5.55]);  % Why there are spikes in range cells where there is no target ??
        title(['rcell = ',num2str(rcell)],'Fontsize',16);

        figure;
        plot(doppler_vel_symm,abs(dopplerfft_matrix(rcell,:)));
        %ylim([0,max_val+10]);
        ylim([0,200.55]);  % Why there are spikes in range cells where there is no target ??
        title(['rcell = ',num2str(rcell)],'Fontsize',16);

        figure;
        plot(doppler_vel_symm,abs(dopplerfft_matrix(rcell,:)));
        %ylim([0,max_val+10]);
        ylim([0,0.15]);  % Why there are spikes in range cells where there is no target ??
        title(['rcell = ',num2str(rcell)],'Fontsize',16);
        pause;
        %
    end
end


dopplerfft_matrix_abs = abs(dopplerfft_matrix);

max_doppler_abs = max(max(dopplerfft_matrix_abs));
min_doppler_abs = min(min(dopplerfft_matrix_abs));
fprintf("max_doppler_abs=%d  min_doppler_abs=%d \n",max_doppler_abs,min_doppler_abs);
dopplerfft_matrix_abs_db = db(dopplerfft_matrix_abs./max(max(dopplerfft_matrix_abs)));


pause;
close all;

% Mesh Grid - doppler vel vs Range - in DB
figure;
[X,Y] = meshgrid(doppler_vel_symm,range_vect(1:RFFT_BY_2));
s = mesh(X,Y,dopplerfft_matrix_abs_db);
xlabel ('Velocity (m/s) ');
ylabel('Range(m)');
s.FaceColor = 'flat';
colorbar;
ylim ([0 RFFT_BY_2]);
xlim([doppler_vel_symm(1),doppler_vel_symm(DFFT)]);
zlim([min(min(dopplerfft_matrix_abs_db)),max(max(dopplerfft_matrix_abs_db))
]);


% Mesh Grid - doppler vel vs Range
figure;
[X,Y] = meshgrid(doppler_vel_symm,range_vect(1:RFFT_BY_2));
s = mesh(X,Y,dopplerfft_matrix_abs);
xlabel ('Velocity (m/s) ');
ylabel('Range(m)');
s.FaceColor = 'flat';
colorbar;
ylim ([0 RFFT_BY_2]);
xlim([doppler_vel_symm(1),doppler_vel_symm(DFFT)]);
zlim([0.3 2e3]);


% Mesh Grid  around zero vel - doppler vel vs Range
figure;
[X,Y] = meshgrid(doppler_vel_symm,range_vect(1:RFFT_BY_2));
s = mesh(X,Y,dopplerfft_matrix_abs);
xlabel ('Velocity (m/s) ');
ylabel('Range(m)');
s.FaceColor = 'flat';
colorbar;
%ylim ([0 RFFT_BY_2]);
ylim([0 10]);
%xlim([doppler_vel_symm(1),doppler_vel_symm(DFFT)]);
xlim([-5 5]);
%zlim([0.5 2e3]);
zlim([5 4e3]);

%pause;
%close all;

figure;
imagesc(doppler_vel_symm,range_vect(1:RFFT_BY_2),(dopplerfft_matrix_abs));
colorbar;
ylabel('Range(m)');
xlabel ('Velocity (m/s) ')
str = {'Range - Velocity - Doppler Spectrum';...
            ['Beam Number = ',num2str(beam_num-1)]}; 
title(str,'FontSize',16);
clim ([-10 10])
%caxis ([1e-8, 4e4]);
ylim ([0 RFFT_BY_2]);
return;


function [bm,ak] = define_window_based_fir_filter(fs) 
    %
    cf1 = 5e3;
    cf2 = 8e6;
    nf = fs/2;
    ncf1 = cf1/nf;
    ncf2 = cf2/nf;
    ncf = [ncf1, ncf2];
    RFFT = 4096;
    order = RFFT/16;
    win = window(@hamming, order);
    bm = fir1(order-1,ncf,'bandpass',win);
    ak = 1;
    %

    %{
    cf2 = 8e6; %aaa
    nf = fs/2;
    ncf2 = cf2/nf;
    order = 256;
    win = window(@hamming, order);
    bm = fir1(order-1,ncf2,'low',win);
    ak=1;
    %}
 end



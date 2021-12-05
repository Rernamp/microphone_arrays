clear;close all;


K = 9; 
c = 343;

fs = 8000;

f = 2000;
phi_sig = 0;
teta_sig = 0;
phi_noise = 50;
teta_noise = 0;
c = 343;
d = c/(2*f);
d = 0.04;
L = 32;
N_f = 50; 
mu = 0.01;
c_noise =35;
p_loc = gen_place_el(sqrt(K),sqrt(K),d,d,1)';

[noise,fs] = audioread('nois.wav');
[signal,fs] = audioread('sound.wav');
osh_in = mean(signal.^2)/mean(noise.^2);
t = 2;
time = 0:1/fs:t-1/fs;
signal = signal(1:t*fs);
noise = noise(1:t*fs);

signal_shift = shift_plane(signal,phi_sig,teta_sig,p_loc,fs);
noise_shift = shift_plane(noise,phi_noise,teta_noise,p_loc,fs);


mu = 1;

time = 0:1/fs:(length(noise)-1)/fs;



phi = -90:90;


for l = 1:length(phi)
    signal_shift = shift_plane(signal,phi(l),teta_sig,p_loc,fs);
    %noise_shift = shift_plane(noise,phi(l),teta_noise,p_loc,fs);
    [y_noise,W_n] = spat_filt_wb_time_lc_lms(noise_shift, L, K, mu);
    [y_sig,W_s] = spat_filt_wb_time_lc_lms(signal_shift, L, K, mu);
    osh_phi(l) = mean(y_sig.^2)/mean(y_noise.^2);
    
   
end
%%
teta = 0:90;
for l = 1:length(teta)
    signal_shift = shift_plane(signal,phi_sig,teta(l),p_loc,fs);
    %noise_shift = shift_plane(noise,phi_noise,teta(l),p_loc,fs);
    [y_noise,W_n] = spat_filt_wb_time_lc_lms(noise_shift, L, K, mu);
    [y_sig,W_s] = spat_filt_wb_time_lc_lms(signal_shift, L, K, mu);
    osh_teta(l) = mean(y_sig.^2)/mean(y_noise.^2);
   
end

%%

SNR_phi = 20*log10(osh_phi-osh_in);
figure()
plot(phi,SNR_phi)
ylabel('ОСШ')
xlabel('Порядок фильтра, K = 4')
title('ОСШ взависимости от порядка фильтра для LC RLS')
grid on
%%

SNR_teta = 20*log10(osh_teta-osh_in);
figure()
plot(teta,SNR_teta)
ylabel('ОСШ')
xlabel('Число микрофонов, J = 32')
title('ОСШ взависимости от числа элементов МР для LC RLS')
grid on

save('SNR_LMS.mat','phi','teta','SNR_phi','SNR_teta');

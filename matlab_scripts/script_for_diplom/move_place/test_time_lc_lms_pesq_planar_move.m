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
p_loc(3,6) = p_loc(3,6) + d/10; 
p_loc(2,6) = p_loc(2,6) - d/20; 
p_cil = p_loc;
[noise,fs] = audioread('nois.wav');
[signal,fs] = audioread('speech_dft_8kHz.wav');
osh_in = mean(signal.^2)/mean(noise.^2);
t = 2;
time = 0:1/fs:t-1/fs;
signal = signal(1:t*fs);
noise = noise(1:t*fs);

signal_shift = shift_plane(signal,phi_sig,teta_sig,p_loc,fs);
noise_shift = shift_plane(noise,phi_noise,teta_noise,p_loc,fs);

sig_in_MR = signal_shift + noise_shift;
mu = 1;

time = 0:1/fs:(length(noise)-1)/fs;

mu = 1;

time = 0:1/fs:(length(noise)-1)/fs;
%%
L_max = 56;
L_min = 16;
L = L_min:L_max; 
i =1; 

for l = L_min:4:L_max
    [y_cl,y,W_n] = spat_filt_wb_time_lc_lms_pesq(signal_shift,sig_in_MR,l,K,mu);
    k_l(i,1:2) = pesqbin(y_cl,y,fs,'nb');
    i = i+1;
   
end
%%

mu = 0.1; 

k = [4 9 16 25];
L_k = 32;
for p = 1:length(k)
    
    p_loc = gen_place_el(sqrt(k(p)),sqrt(k(p)),d,d,1)'; 
    signal_shift = shift_plane(signal,phi_sig,teta_sig,p_loc,fs);
    noise_shift = shift_plane(noise,phi_noise,teta_noise,p_loc,fs);
    
    sig_in_MR = signal_shift + noise_shift;    
    [y_cl,y,W_n] = spat_filt_wb_time_lc_lms_pesq(signal_shift,sig_in_MR,L_k,k(p),mu);
    k_K(p,1:2) = pesqbin(y_cl,y,fs,'nb');
    
end
%%
l = L_min:4:L_max;

figure()

plot(l,k_l(:,1))
grid on
title("PESQ для LC RLS")
xlabel("J , K = 4")
ylabel("PESQ")



figure()

plot(k,k_K(:,1))
grid on
title("PESQ для LC RLS")
xlabel("K, J = 32")
ylabel("PESQ")

save('PESQ_LMS_move.mat','l','K','k_l','k_K');

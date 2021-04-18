clear;close all;

c = 343;
f_re = 4000;
lamda = c/f_re;
K = 4;
f_up = 4000;   
f_down = 0;  
f_den = 0.5;   
L = 30;
d = 0.04;

thetaad = 0;     
thetaan = 60;      

[noise,fs] = audioread('nois.wav');

[sig,fs] = audioread('speech_dft_8kHz.wav');

t = 2;
%%
sig = sig(1:t*fs);
noise = noise(1:t*fs);

sig_in_MR_noise = zeros(K, length(sig));
sig_in_MR_sig = zeros(K, length(sig));

tau = (d*sind(thetaan))/(c);
osh_in = mean(sig.^2)/mean(noise.^2);
%%
for k = 1:K
    [y_shift] = shift( noise, tau*(k-1), fs);
    sig_in_MR_noise(k,:) = y_shift + sig';
    sig_in_MR_sig(k,:) = sig;
end

mu = 1;

time = 0:1/fs:(length(noise)-1)/fs;
%%
L_max = 64;
L_min = 16;
L = L_min:L_max; 
 

for l = L_min:2:L_max
    [y_cl,y,W_n] = spat_filt_wb_time_lc_lms_pesq(sig_in_MR_sig,sig_in_MR_noise,l,K,mu);
    k_l(l-L_min+1,1:2) = pesqbin(y_cl,y,fs,'nb');
   
end
%%

mu = 0.1; 

k_max = 8;
k_min = 2;
L_k = 32;
for p = k_min:k_max
    sig_in_MR_noise = zeros(p, length(sig));
    sig_in_MR_sig = zeros(p, length(sig));
    for n = 1:p
        [y_shift] = shift( noise, tau*(n-1), fs);
        sig_in_MR_noise(n,:) = y_shift+sig';
        sig_in_MR_sig(n,:) = sig;
        
    end
    [y_cl,y,W_n] = spat_filt_wb_time_lc_lms_pesq(sig_in_MR_sig,sig_in_MR_noise,L_k,p,mu);
    k_K(p-k_min+1,1:2) = pesqbin(y_cl,y,fs,'nb');
    
end
%%
l = L_min:L_max;

figure()

plot(l,k_l(:,1))
grid on
title("PESQ")
xlabel("J, K = 4")
ylabel("ОСШ")

p = k_min:k_max;

figure()

plot(p,k_K(:,1))
grid on
title("PESQ")
xlabel("K, J = 32")
ylabel("ОСШ")
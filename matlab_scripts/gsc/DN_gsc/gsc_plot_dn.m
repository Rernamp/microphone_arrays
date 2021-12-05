clear;
close all;

N = 9; 
c = 343;

fs = 8000;

f = 2000;
phi_sig = 0;
teta_sig = 0;
phi_noise = 60;
teta_noise = 0;
c = 343;
J = 16;
% d = c/(2*f);
t = 2;
time = 0:1/fs:t-1/fs;
d = 0.04;
N_f = 50; 
mu = 0.1;
del_p = 0.005;
p_loc = gen_place_el(sqrt(N),sqrt(N),d,d,1)';
% p_loc(2,3) = p_loc(2,3) + del_p;
% p_loc(3,3) = p_loc(3,3) - del_p;
% % 
% p_loc(2,6) = p_loc(2,6) - del_p;
% p_loc(3,6) = p_loc(3,6) + del_p;
[noise,fs] = audioread('nois.wav');
[signal,fs] = audioread('sound.wav');

t = 2;

signal = signal(1:t*fs);
noise = noise(1:t*fs);

signal_shift = shift_plane(signal,phi_sig,teta_sig,p_loc,fs);
noise_shift = shift_plane(noise,phi_noise,teta_noise,p_loc,fs);

sig_in_MR = signal_shift+noise_shift;
%% algoritm
[y_GSC,w] = GSC_RLS(sig_in_MR,J,N);

B = eye(N-1,N) - [zeros(N-1,1) eye(N-1,N-1)];
f = [1; zeros(J-1,1)];
C = zeros(N*J,J);

for j = 1:J
    C(:,j) = [zeros(1,(j-1)*N) ones(1,N) zeros(1,J*N-j*N)].';
end
     %%
w_q = C*inv(C'*C)*f;

%%
W_fft = reshape(w,N-1,J);
B_B = B'* inv(B*B');
W = B_B*W_fft;
W_w = reshape(W,N*J,1);
W_W = w_q - W_w ;
%%
[BP_PHO,BP_TETA] = plot_bp_for_place(p_loc,phi_sig,teta_sig,50,W_W,fs,J,N);
%%
figure()
hold on
grid on
plot(time,sig_in_MR(1,:))

% plot(time,y_LMS)
plot(time,y_GSC)
legend('input','LC NLMS','GSC LMS')
% legend('input','GSC LMS')

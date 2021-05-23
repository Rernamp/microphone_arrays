clear;
close all;

N = 16; 
c = 343;
L = 32;
fs = 8000;
f = 2000;
phi_sig = 0;
teta_sig = 0;
phi_noise = 30;
teta_noise = 0;
c = 343;
d = c/(2*f);
d = 0.04;
p_loc = gen_place_el(4,4,d,d,1)';

[noise,fs] = audioread('nois.wav');

[sig,fs] = audioread('speech_dft_8kHz.wav');

t = 2;

time = 0:1/fs:(length(noise)-1)/fs;
%%
sig = sig(1:t*fs);
noise = noise(1:t*fs);
power_sig = mean(sig.^2);
power_noise = mean(noise.^2);

mu = 1;



signal_shift = shift_plane(sig,phi_sig,teta_sig,p_loc,fs);
noise_shift = shift_plane(noise,phi_noise,teta_noise,p_loc,fs);

sig_noise = signal_shift + noise_shift;

[y,W] = spat_filt_wb_time_lc_lms(sig_noise, L, N, mu);
 
%%
phi_const = 0;
teta_const = 30;

T_S = 1/fs;

f = [0:L/2-1 -L/2:-1]*fs/L;   
f_pos = f( 1 : length(f)/2 ); 
phi = -90:0.5:90;
teta = 0:0.5:90;    
a_phi = [-cosd(teta).*cosd(phi_const) ; -cosd(teta).*sind(phi_const) ; -sind(teta)];
a_teta = [-cosd(teta_const).*cosd(phi) ; -cosd(teta_const).*sind(phi) ; -sind(teta_const).*ones(1,length(phi))];
    

   
W_fft = reshape(W,N,L);
% teta const
for f_i = 1:length(f_pos)
    for i_phi = 1:length(phi)
        a_teta = [-cosd(teta_const).*cosd(phi(i_phi)) ; -cosd(teta_const).*sind(phi(i_phi)) ; -sind(teta_const)];
        tau = a_teta'*p_loc/c;
        d = exp(-1i*f_pos(f_i)*(tau))';
        for k =2:L
            d = [d ; exp(-1i*f_pos(f_i)*(tau + (k-1)*T_S))'];
        end
        BP_teta(f_i,i_phi) = W'*d;
    %     d = ;

    end
end
%%
%phi const
for f_i = 1:length(f_pos)
    for i_teta = 1:length(teta)
        a_phi = [-cosd(teta(i_teta)).*cosd(phi_const) ; -cosd(teta(i_teta)).*sind(phi_const) ; -sind(teta(i_teta))];
        tau = a_phi'*p_loc/c;
        d = exp(-1i*f_pos(f_i)*(tau))';
        for k =2:L
            d = [d ; exp(-1i*f_pos(f_i)*(tau + (k-1)*T_S))'];
        end
        BP_phi(f_i,i_teta) = W'*d;
    %     d = ;

    end
end

BP_phi = abs(BP_phi).^2; 
BP_teta = abs(BP_teta).^2;
    
BP_phi_dB = 10*log10(BP_phi);
BP_teta_dB = 10*log10(BP_teta);
    
BP_phi_dB(BP_phi_dB < -40) = -40;
BP_teta_dB(BP_teta_dB < -40) = -40;


    [X_phi,Y_phi] =  meshgrid(f_pos,teta);
    
    figure()
    tit = strcat("Const \phi = " ,string(phi_const));
    surf(X_phi,Y_phi,BP_phi_dB')
    grid on
    zlim([-40 50])
    xlabel("frequency , f");
    title(tit)
    ylabel("angle,\theta");
    zlabel("BP, dB");
    shading interp 
    grid on %
%     colormap gray
    
    [X_teta,Y_teta] = meshgrid(f_pos,phi);
    
    figure()
    tit = strcat("Const \theta = " ,string(teta_const));
    surf(X_teta,Y_teta,BP_teta_dB')
    grid on
    zlim([-40 50])
    xlabel("frequency , f");
    ylabel("angle,\phi");
    title(tit)
    zlabel("BP, dB");
    shading interp 
    grid on %
    %%
    [B,BB] = plot_bp_for_place(p_loc,phi_const,teta_const,L/2,W,fs,L,N);
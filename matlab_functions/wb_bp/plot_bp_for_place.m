function [BP_phi_dB,BP_teta_dB] = plot_bp_for_place(p_el,phi_const,teta_const,N_f,W,fs,J,N_el)
        
    c = 343;
    fs = 8e3;
    f = [0:N_f/2-1 -N_f/2:-1]*fs/N_f;   
    f_pos = f( 1 : length(f)/2 ); 

    phi = -90:0.5:90;
    teta = 0:0.5:90;
    W_fft = reshape(W,N_el,J);
    W_fft = fft(W_fft,N_f,2);
   
    
%      a_phi = [-sind(teta).*cosd(phi_const) ; -sind(teta).*sind(phi_const) ; -cosd(teta)];
%      a_teta = [-sind(teta_const).*cosd(phi) ; -sind(teta_const).*sind(phi) ; -cosd(teta_const).*ones(1,length(phi))];

     a_phi = [-cosd(teta).*cosd(phi_const) ; -cosd(teta).*sind(phi_const) ; -sind(teta)];
     a_teta = [-cosd(teta_const).*cosd(phi) ; -cosd(teta_const).*sind(phi) ; -sind(teta_const).*ones(1,length(phi))];
    
    a_p_phi = a_phi'*p_el;
    a_p_teta = a_teta'*p_el;
    
    for f_i = 1:length(f_pos)
%         BP_phi(f_i,:) = (W_fft(:,f_i)')*exp(-1i*2*pi*f_pos(f_i).*a_p_phi'/c);
%         BP_teta(f_i,:) = (W_fft(:,f_i)')*exp(-1i*2*pi*f_pos(f_i).*a_p_teta'/c);
        BP_phi(f_i,:) = exp(-1i*2*pi*f_pos(f_i).*a_p_phi/c)*(W_fft(:,f_i));
        BP_teta(f_i,:) = exp(-1i*2*pi*f_pos(f_i).*a_p_teta/c)*(W_fft(:,f_i));
    end   
    BP_phi = abs(BP_phi).^2; 
    BP_teta = abs(BP_teta).^2;
    
    BP_phi_dB = 10*log10(BP_phi);
    BP_teta_dB = 10*log10(BP_teta);
    
    BP_phi_dB(BP_phi_dB < -40) = -40;
    BP_teta_dB(BP_teta_dB < -40) = -40;
    
    [X_phi,Y_phi] =  meshgrid(teta,f_pos);
    
    figure()
    tit = strcat("Фиксированый угол \phi = " ,string(phi_const));
    surf(X_phi,Y_phi,BP_phi_dB)
    grid on
 
    ylabel("Частота , Гц");
    %title(tit)
    xlabel("Угол подъёма, град");
    zlabel("BP, дБ");
    shading interp 
    grid on %
    colormap gray
    
    [X_teta,Y_teta] = meshgrid(phi,f_pos);
    
    figure()
    tit = strcat("Фиксированый угол \theta = " ,string(teta_const));
    surf(X_teta,Y_teta,BP_teta_dB)
    grid on
  
    ylabel("Частота , Гц");
    xlabel("Угол азимута, град");
    %title(tit)
    zlabel("BP, дБ");
    shading interp 
    grid on %
    colormap gray
%     
end


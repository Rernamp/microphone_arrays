function [BP_phi_dB,BP_teta_dB] = wb_bp_without_fft(p_loc,phi_const,teta_const,N_f,W,fs,L)

        
    c = 343;
    fs = 8e3;
    f = [0:N_f/2-1 -N_f/2:-1]*fs/N_f;   
    f_pos = f( 1 : length(f)/2 ); 
    T_S = 1/fs;

    phi = -90:0.5:90;
    teta = 0:0.5:90;
    
    for f_i = 1:length(f_pos)        
        for i_phi = 1:length(phi)
            a_teta = [-cosd(teta_const).*cosd(phi(i_phi)) ; -cosd(teta_const).*sind(phi(i_phi)) ; -sind(teta_const)];
            tau = a_teta'*p_loc/c;
            d = exp(-1i*2*pi*f_pos(f_i)*(tau));
            for k_i = 2:L

                d = [d exp(-1i*2*pi*f_pos(f_i)*(tau + (k_i-1)*T_S))];
            end
            BP_teta(f_i,i_phi) = d*W;
        end
    end

    %%

    for f_i = 1:length(f_pos)        
        for i_teta = 1:length(teta)
                a_phi = [-cosd(teta(i_teta)).*cosd(phi_const) ; -cosd(teta(i_teta)).*sind(phi_const) ; -sind(teta(i_teta))];
                tau = a_phi'*p_loc/c;

                d = exp(-1i*2*pi*f_pos(f_i)*(tau));
                for k_i = 2:L

                    d = [d exp(-1i*2*pi*f_pos(f_i)*(tau + (k_i-1)*T_S))];
                end
                BP_phi(f_i,i_teta) = d*W;

        end
    end
    %%
    BP_phi = abs(BP_phi).^2; 
    BP_teta = abs(BP_teta).^2;

    BP_phi_dB = 10*log10(BP_phi);
    BP_teta_dB = 10*log10(BP_teta);

    BP_phi_dB(BP_phi_dB < -40) = -40;
    BP_teta_dB(BP_teta_dB < -40) = -40;


    [X_phi,Y_phi] =  meshgrid(teta,f_pos);

    figure()
    tit = strcat("Const \phi = " ,string(phi_const));
    surf(X_phi,Y_phi,BP_phi_dB)
    grid on
    zlim([-40 50])
    ylabel("frequency , f");
    title(tit)
    xlabel("angle,\theta");
    zlabel("BP, dB");
    shading interp 
    grid on %
    %     colormap gray

    [X_teta,Y_teta] = meshgrid(phi,f_pos);

    figure()
    tit = strcat("Const \theta = " ,string(teta_const));
    surf(X_teta,Y_teta,BP_teta_dB)
    grid on
    zlim([-40 50])
    ylabel("frequency , f");
    xlabel("angle,\phi");
    title(tit)
    zlabel("BP, dB");
    shading interp 
    grid on %

end
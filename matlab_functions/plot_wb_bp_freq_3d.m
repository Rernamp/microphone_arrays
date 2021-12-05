function [wfreq_h] = plot_wb_bp_freq_3d(W, Nel, J, pos, num_bin, azi_plot, elv_plot )
    c = 343;
    fs = 8e3;

    f = [0:num_bin/2-1 -num_bin/2:-1]*fs/num_bin;   
    f_pos = f( 1 : length(f)/2  )   
    azi = -90:0.5:90;
    elv = 0:0.5:90;
    P_azi =  zeros( length(azi)  , 1 );
    BP_azi = zeros( length(f_pos), length(azi) );
    P_elv =  zeros( length(elv)  , 1 );
    BP_elv = zeros( length(f_pos), length(elv) );
    tic
    wfreq_h = fft(reshape(W,Nel,J), num_bin, 2);  
    for ff = 1 : length(f_pos)
        k_num = 2*pi*f_pos(ff)/c;
        for aa = 1 : length(azi)
            k = -k_num*[ cosd(elv_plot)*cosd(azi(aa));
                         cosd(elv_plot)*sind(azi(aa));
                         sind(elv_plot)];               
            P_azi(aa) =  exp(-1i*k.'*pos)* wfreq_h(:, ff);                
        end
        for ee = 1 : length(elv)
            k = -k_num*[ cosd(elv(ee))*cosd(azi_plot);
                         cosd(elv(ee))*sind(azi_plot);
                         sind(elv(ee))];               
            P_elv(ee) =  exp(-1i*k.'*pos)* wfreq_h(:, ff);                
        end
        BP_azi(ff,:) = abs(P_azi).^2;
        BP_elv(ff,:) = abs(P_elv).^2;
    end

    BP_azi_dB = 10*log10(BP_azi);
    BP_elv_dB = 10*log10(BP_elv);

    BP_azi_dB(BP_azi_dB < -40) = -40;
    BP_elv_dB(BP_elv_dB < -40) = -40;
    toc



     %% 
    [AZI,FRREQ_AZI] = meshgrid(azi, f_pos);
%     BP_azi_dB = floor(BP_azi_dB*1e3)/1e3;
    figure1 = figure();
    
    surf(AZI,FRREQ_AZI,BP_azi_dB); 
    
    %%
    [ELV,FRREQ_ELV] = meshgrid(elv, f_pos);
%     BP_elv_dB = floor(BP_elv_dB*1e1)/1e1;
    figure2 = figure();
   
    surf(ELV,FRREQ_ELV,BP_elv_dB); 
   
end


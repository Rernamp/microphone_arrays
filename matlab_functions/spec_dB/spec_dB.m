function [  ] = spec_dB( x_in , fs)

    x_fft = fft(x_in);
    
    ahc_x = mag2db(abs(x_fft));
    
    f = (0:length(x_in)-1)*fs/length(x_in);
    
    plot(f, ahc_x);
    set(0,'DefaultAxesFontSize',14,'DefaultAxesFontName','Times New Roman');
    set(0,'DefaultTextFontSize',14,'DefaultTextFontName','Times New Roman'); 
    grid on;
    xlabel('„астота, Hz');
    ylabel("BP, дЅ")
end


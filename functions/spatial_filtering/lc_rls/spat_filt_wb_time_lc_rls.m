function [ y,W ] = spat_filt_wb_time_lc_rls( s_in, J, M)

    delta_min_2 = 1000;
lamda = 1;
%инициализация
x = zeros(M*J,1);
R_inv = delta_min_2*eye(M*J,M*J);

f = [1; zeros(J-1,1)];

C = [ones(M,1) ;zeros(M*(J-1),1)];
for l = 1 : J-1
    C = [C [zeros(l*M,1) ; ones(M,1) ; zeros(M*(J-l-1),1)]];
end

Gama = R_inv * C; 

Fi_inv = inv(C' * Gama);

W = Gama * Fi_inv * f;
y = zeros(1,length(s_in(1,:)));
for k = 1:length(s_in(1,:))
    
    x = [s_in(:,k); x(1:end-M)];
    y(k) = W' * x;
    g = (R_inv * x)./(lamda + x' * R_inv * x);
    R_inv = (R_inv - g * x' * R_inv)./(lamda);
    
    nu = C' * g;
 
    n = x' * Gama;
  
    
    om = (Fi_inv * nu)./(1 - n * Fi_inv * nu);
    
    Fi_inv = lamda.*(Fi_inv + om * n * Fi_inv);
    
    Gama = (Gama - g * n)./(lamda);
    
    alpha = 0 - W' * x;
    
    W = W + (g - lamda * Gama * om).* conj(alpha);
    
end

end


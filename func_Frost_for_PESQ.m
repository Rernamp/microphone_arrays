function [y_clear,y,W] = func_Frost_for_PESQ(s_in_clear,s_in, J, K, mu)
    %создаю начальные и вспомогательные вектора
    y = zeros(length(s_in(1,:)),1);
    y_clear = zeros(length(s_in(1,:)),1);
    f = [1 ;zeros(J-1,1)];
    C = [ones(K,1) ;zeros(K*(J-1),1)];
    for l = 1 : J-1
        C = [C [zeros(l*K,1) ; ones(K,1) ; zeros(K*(J-l-1),1)]];
    end
      
    W_p = C*inv(C'*C)*f;       
    
    W = zeros(K, J);
    %в статье у них получаетс€ вектор размером KJ на 1. ѕоэтому разбиваю
    %дл€ себ€.  онечно можно оставить и в таком варианте,но нужно изменить
    %код
    for b = 1:K
        W(b,:) = W_p(1+(b-1)*J:J*b);
    end
    %вектора дл€ посчета суммы
    I_for_W = ones(1,K);
    I_for_x = ones(K,1);
    %цикл фильтрации
    x_i = zeros(J, K);
    x_i_clear = zeros(J, K);
    for k = 1 : length(s_in)
        
        %набор новых значений
        x_i = [s_in(:, k).' ; x_i(1:J - 1, :)];
        x_i_clear = [s_in_clear(:, k).' ; x_i_clear(1:J - 1, :)];
        %цикл фильтрации
        for z = 1 : K
            y(k) = y(k) + (W(z,:) * x_i(:,z));
            y_clear(k) = y_clear(k) + (W(z,:) * x_i_clear(:,z));
        end
        
        %цикл дл€ алгоритма
        
        Dop = I_for_W * W - mu * y(k) * I_for_x' * x_i';
        for n = 1 : J
                        
            for l = 1:K
                W (l, n) = W(l, n) - mu * y(k) .* x_i(n, l).' - (Dop(n) /K) + f(n)/K;
            end
                       
        end
        
    end
end


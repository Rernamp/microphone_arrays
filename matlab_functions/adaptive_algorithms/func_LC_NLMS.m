function [ y,W ] = func_LC_NLMS(s_in, J, M, mu )

x = zeros(M*J,1);
f = [1; zeros(J-1,1)];

delta_min_2 = 2.2204460492503131e-016;

C = [ones(M,1) ;zeros(M*(J-1),1)];
for l = 1 : J-1
    C = [C [zeros(l*M,1) ; ones(M,1) ; zeros(M*(J-l-1),1)]];
end
    
Q = C * inv(C'*C);

h = Q*f;
y = zeros(1,length(s_in));
for i = 1:length(s_in) 

    x = [s_in(:,i); x(1:end-M)];
    y(i) = h' * x;
    
    alpha = - h' * x;
    k_k = ((x'*x + (x'*Q)*(C'*x) + delta_min_2));
    alpha = alpha*inv(k_k);
    h = h + 2*mu*x*alpha;
    
    h = h + Q*(f-C'*h);
    
end
W = h;

end


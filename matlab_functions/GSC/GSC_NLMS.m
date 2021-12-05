function [ y ] = GSC_NLMS(s_in, J, N, mu)
    
    f = [1; zeros(J-1,1)];
    B = eye(N-1,N) - [zeros(N-1,1) eye(N-1,N-1)];

    C = zeros(N*J,J);

     for j = 1:J
        C(:,j) = [zeros(1,(j-1)*N) ones(1,N) zeros(1,J*N-j*N)].';
     end
     %%
     w_q = C*inv(C'*C)*f;
     x = zeros((N-1)*J,1);
%      h_a = w_q(1:(N-1)*J);
     h_a = zeros((N-1)*J,1);
     x_w_q = zeros(N*J,1);
     delta_min_2 = 2.2204460492503131e-016;
     y = zeros(1,length(s_in(1,:)));
     d = zeros(1,length(s_in(1,:)));

     for i = 1:length(s_in) 
        x_B = B*s_in(:,i);
        x_w_q = [s_in(:,i); x_w_q(1:end-N)];
        x = [x_B; x(1:end-N+1)];
        y(i) = h_a' * x;

        d(i) = w_q'*x_w_q;
        alpha = d(i) - y(i);
        k_k = ((x'*x + delta_min_2));
        alpha = alpha*inv(k_k);
        h_a = h_a + 2*mu*x*alpha;

        y(i) = d(i) - y(i);
     end
end
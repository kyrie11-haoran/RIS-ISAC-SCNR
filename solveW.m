function [W] = solveW(hB,hR,G,P_0,R,N_t,C,K,L,lambda_c,lambda_s,f_k,f_K1,v,sigma_k2,sigma_s2,a_c,a_s,alpha_c,alpha_s,epsilon)
   
%%%%%  参数设置
  
   h_tilde = zeros(N_t+1,K);%%  预分配
   A = zeros(N_t+1,N_t+1,K);
   B = zeros(N_t+1,N_t+1,K);
  a_c_tilde = zeros(N_t+1,C);
   E = zeros(N_t+1,N_t+1,C);

   Theta = diag(v);
   h_k = zeros(N_t, K); % 初始化 h_k 
for k=1:K
    h_k(:,k) = (hB(:,k)' + (hR(:,k)' * Theta * G))';    %% theta is random
end
   for i = 1:K
        % 定义 A_k (24)
        h_tilde(:,i) = [h_k(:,i)',0].';  
        A(:,:,i) = h_tilde(:,i) * h_tilde(:,i)';    

        % 定义 B_k (23)
        B(:,:,i) = [abs(f_k(i))^2 * h_k(:,i) * h_k(:,i)',  zeros(N_t,1);  
                    2 * f_k(i).' * h_k(:,i)', 0];
   end
    for c = 1:C
       %定义 E_c
      a_c_tilde(:,c) = [a_c(:,c)', 0].'; 
        E(:,:,c) = (alpha_c)^2 * a_c_tilde(:,c) * a_c_tilde(:,c)';
    end
    
     %定义 D
       D = [abs(f_K1)^2 * (alpha_s)^2 * a_s * (a_s)',  zeros(N_t,1);  
           2 * conj(f_K1) * (a_s)', 0];

       
 cvx_begin sdp 
 variable W_jian(N_t+1,N_t+1,K) hermitian semidefinite
 variable c_k(1,K)  nonnegative
 variable c_s  nonnegative


    %%主函数
   %minimize(square_pos(norm(sum(W_jian(1:N_t, 1:N_t), 3) - R, 'fro')))
     maximize(lambda_c * (sum(c_k)) + lambda_s * c_s)%%(34a)
    subject to
    %%限制条件(34b)
    %%(21)
for i = 1:K
    real(trace(B(:,:,i) * W_jian(:,:,i)) - abs(f_k(i))^2 * (sigma_k2 + trace(A(:,:,i) * sum(W_jian,3))))  >=  c_k(i);%+ 1/6 * log(2)^3 * pow_pos(c_k(i), 3) + 1/24 * log(2)^4 * pow_pos(c_k(i), 4) + 1/120 * log(2)^5 * pow_pos(c_k(i), 5);   
end
    %%(22)
   real(trace(D * sum(W_jian,3)) - abs(f_K1)^2 * (sigma_s2 + trace(sum(E,3) * sum(W_jian,3)))) >= c_s; %+ 1/6 * log(2)^3 * pow_pos(c_s, 3) + 1/24 * log(2)^4 * pow_pos(c_s, 4) + 1/120 * log(2)^5 * pow_pos(c_s, 5);%%

    %%(28)(29)(30)
   square_pos(norm(sum(W_jian(1:N_t, 1:N_t), 3) - R, 'fro')) <= epsilon;          %(28)
 
    %diag(sum(W_jian(1:N_t, 1:N_t),3)) == (P_0 * ones(N_t , 1))/N_t; 
   
    trace(sum(W_jian(1:N_t, 1:N_t),3))<=P_0;
    %(29)
%(lambda_c * (sum(c_k)) + lambda_s * c_s)<=20;
    cvx_end
    beampattern_constraint = norm(sum(W_jian(1:N_t, 1:N_t), 3) - R, 'fro')^2;
    disp(beampattern_constraint);

    W = W_jian;

function [v]  = solvev(N_t,L,G,hR,hB,w,K,f_k,sigma_k2,lambda_c)%%%%% N_t,L,G,hR,hB,w,K,f_k,sigma_k2,lambda_c
iter = 1;
iter_max = 20;
tau = 1.5;  %下一次惩罚的倍数
tau_1 = 0.8;%v-v_before的门限
tau_2 = 1.2;%xi的模
rho_max = 100;%最大惩罚倍
rho = 10;%初始惩罚倍数
data = load('v_storage.mat'); 
v_before =data.v_storage(:, 2);


while(iter <= iter_max)
 cvx_begin 
 variable v(L,1) complex
 variable xi(2*L,1) nonnegative
 variable c_k(1,K)  nonnegative
 variable p(1,K)
 Sigma = zeros(N_t,N_t,K);
for k = 1:K
    for i = 1:K
        if i ~= k
    Sigma(:,:,k) = Sigma(:,:,k)+ w(:,i) * w(:,i)';
        end
    end
end

    Q = zeros(L, K);
    M = zeros(L, L, K);
    z = zeros(1,K);


for k = 1:K    
    
    Q(:,k) = abs(f_k(k))^2 * diag(hR(:,k)) * G * Sigma(:,:,k) * conj(hB(:,k)) - f_k(k)' * diag(hR(:,k)) * G * w(:,k); 

    M(:,:,k) = abs(f_k(k))^2 * diag(hR(:,k)) * G * Sigma(:,:,k) * G' * diag(hR(:,k))';
    
    z(k)  =  2 * real(conj(f_k(k)) * (hB(:, k).' * w(:, k))) ...
             - abs(f_k(k))^2 * real((sigma_k2 + hB(:, k).' * Sigma(:,:,k) * conj(hB(:, k))));   
end
     

%%主函数

maximize(lambda_c*sum(c_k)- rho * norm(xi,1))
%maximize ( sum(t) - rho * norm(xi,1))   
%%约束条件    
    subject to
%%(41)
for k=1:K
    
    norm([sqrt(M(:,:,k)) * conj(v) ; (1 - p(k) )/2 + real(Q(:,k)' * conj(v))] , 2 ) <= (1 + p(k)  )/2 - real(Q' * conj(v));
  p(k) <= z(k) - c_k(k);

end

%%(31)

for l = 1:L
    conj(v(l)) * v(l) <= 1 + xi(l); 
end
%%(32)
for l = 1:L
  conj(v_before(l)) * v_before(l) - 2 * real(v(l)' .* v_before(l)) <= -1 + xi(l+L);
end

   cvx_end

  rho = min(tau * rho , rho_max);
  iter = iter + 1;
    
if  (norm(v - v_before, 1) <= tau_1) && (norm(xi, 1) <= tau_2)
   disp(norm(v - v_before, 1));
   disp(norm(xi, 1));
    break;% End while loop
end
   disp(['v-v is ',num2str(norm(v - v_before, 1))]);
   disp(['norm(xi) is ',num2str(norm(xi, 1))]);
   disp(v);
   v_before = v;
end
   
end

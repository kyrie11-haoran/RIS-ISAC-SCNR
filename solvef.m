function [f] = solvef(G, w, hB, hR, v, sigma_k2, K)

    f = zeros(K,1);
        for k = 1:K
            denominator3 =0;
            numerator3 = (hB(:,k)' + v.' * diag(hR(:,k)) * G) * w(:, k);
            
            for i = 1:K
                if i ~= k
                    denominator3 = denominator3 + abs((hB(:,k).' + v.' * diag(hR(:,k)) * G) * w(:, i)).^2;
                end
            end
            
            f(k) = numerator3 / denominator3+sigma_k2;
        end

    
end
 





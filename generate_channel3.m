function [hR3, hB3 , G3] = generate_channel3(PL_BR,PL_Rk, PL_Bk, Rician_factor_dB_1,Rician_factor_dB_2, L3, N_t,K)
    % 小尺度衰落 rician信道生成，NOS+LNOS 
    G_LOS =  diag(sqrt( PL_BR * Rician_factor_dB_1/(Rician_factor_dB_1 + 1))) * ones(L3,N_t); % 直路径信号的幅度

    G_NLOS =  diag(sqrt(PL_BR / (Rician_factor_dB_1 + 1))) * sqrt(1/2)*(randn(L3,N_t) + 1i * randn(L3,N_t)) ; % 多径信号的幅度

    G3 = G_NLOS +  G_LOS ;

    % LoS 部分 (RIS 到用户)
    hR_LoS =   ones(L3, K)*diag(sqrt(PL_Rk * Rician_factor_dB_1 / (Rician_factor_dB_1 + 1))); % 直路径信号的幅度

    % NLoS 部分 (RIS 到用户)
    hR_NLoS =    (randn(L3, K) + 1i * randn(L3, K))*diag(sqrt(PL_Rk /2*(Rician_factor_dB_1 + 1))); % 多径信号的幅度

    % 组合 LoS 和 NLoS 成分
    hR3 = hR_LoS + hR_NLoS;

    % NLoS 部分 (BS 到用户) - hB 只有 NLoS 成分
    hB_NLoS =  (randn(N_t,K) + 1i * randn(N_t,K))*diag(sqrt(PL_Bk / 2*(Rician_factor_dB_2 + 1))); % 多径信号的幅度

    hB3 = hB_NLoS;
end
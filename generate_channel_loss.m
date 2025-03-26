function [PL_BR,PL_Rk,PL_Bk,PL_Bt,PL_Rt] = generate_channel_loss(dBk,dBR,dRk,dBt,dRt, delta_0, alpha_Rk,alpha_BR,alpha_Bk,alpha_Bt,alpha_Rt)
% 大尺度衰落 PL(d)=delta_0 * (d/d_0)^(-alpha) ,set d_0=1

PL_Rk = (delta_0) * (dRk).^(-alpha_Rk);
PL_BR = (delta_0) * (dBR).^(-alpha_BR);
PL_Bk = (delta_0) * (dBk).^(-alpha_Bk);
PL_Bt = (delta_0) * (dBt).^(-alpha_Bt);
PL_Rt = (delta_0) * (dRt).^(-alpha_Rt);

end

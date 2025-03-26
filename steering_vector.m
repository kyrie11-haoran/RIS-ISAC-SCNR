function [a_s ] = steering_vector(Phi_s , N_t )
   %默认d =lamda/2
   a_s = zeros(N_t,1);
for s = 1:N_t
    a_s(s) = exp(1j * pi * (s - 1) * sin(Phi_s));
end

end
  
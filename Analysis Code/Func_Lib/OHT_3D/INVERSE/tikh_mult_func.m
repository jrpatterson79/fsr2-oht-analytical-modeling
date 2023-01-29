function [A_prod] = tikh_mult_func(s,transarg,H,C_dd_inv_pf,L)

if transarg == 'notransp'
    A_prod = [C_dd_inv_pf*H*s; L*s];
elseif transarg == 'transp'
    A_prod = (H'*C_dd_inv_pf')*s(1:m,:) + L'*s(((m+1):end),:);
end
end
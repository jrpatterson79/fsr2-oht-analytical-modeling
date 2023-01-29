function [nll] = struct_likelihood(y,X,R,Q,H)

PSI = H*Q*H' + R;
PSI_inv = inv(PSI);
PHI = H*X;
PPP = PHI'*PSI_inv*PHI;
PPP_inv = inv(PPP);

nll = 0.5*sum(log(eig(PSI))) + 0.5*sum(log(eig(PPP))) + ...
    .5*(y'*(PSI_inv - PSI_inv*PHI*(PPP_inv)*PHI'*PSI_inv)*y);



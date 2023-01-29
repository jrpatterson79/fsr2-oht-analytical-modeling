function [obj_func_val] = tikh_obj_eval(y,s,L,C_dd_inv_pf,h_function)

%tikh_obj_eval: Function which evaluates the tikhonov objective function
%given the square root of the inverse of a data covariance matrix (i.e., 1/
%the standard deviations of the data), and a roughening matrix L

%TODO - Make arguments to minres optional
%TODO - Add options for methods for performing the inversions step of Q

resid = y - h_function(s);
obj_resid_vec = C_dd_inv_pf*resid;
obj_reg_vec = L*s;

obj_func_val = norm(obj_resid_vec) + norm(obj_reg_vec);

function [s_best, H_local, objfunc_value] = ql_tikhonov_inv(y,s_init,C_dd,L,h_function,H_tilde_func, varargin)

%ql_tikhonov_inv.m: Inverse problem solver using Tikhonov regularization
%with a generalized roughening matrix L
% Syntax:
% [s_best, H_local, objfunc_value] =
% ql_tikhonov_inv(y,s_init,C_dd,L,forward_function,sensitivity_function,...)
% In everything below:
% m = number of observations
% n = number of parameter values
% p = number of parameters in geostatistical trend (1)
% s_init should be a n by 1 list of initial parameter value guesses
% y should be a  m by 1 list of observed values
% L is thi tikhonov roughening matrix (n x n, preferably sparse)
% C_dd should be an m by m matrix of measurement error
% h_function should be an anonymous function with one argument (input of
% parameter values), and return a vector (value of all observations)
% h_function should be an anonymous function with one argument (input of
% parameter values) and return a matrix (m x n)
%
% Any other variables (i.e., those necessary for evaluation of
% sensitivity_function or forward_function) can be passed after the
% sensitivity function argument as follows:
% ...,'fem',fem,'max_linsearch',10,...
%
% Arguments used by ql_geostat_inv which can be passed in this section:
%    -max_linsearch: Maximum number of steps in linesearch optimization
%    -max_gradevals: Maximum number of gradient evaluations (maximum number
%    of iterations, in essence).
%    -tol_objfunc: Relative tolerance in the objective function change
%    (percent decrease)
%    -tol_s: Relative tolerance in the change in the value of s (root mean
%    square difference from last iteration).
%    -tol_linsearch: Relative tolerance in the change of the objective
%    function during the linesearch optimization (if used)
%

%Create global variables - in case function crashes due to memory or other
%errors, these variables can still be accessed outside of this program.
global H_tilde
global s_hat
start_time = clock;

%Default values for all optimization constants
max_linsearch = 10;
max_gradevals = 20;
tol_objfunc = .001;
tol_s = .001;
tol_linsearch = .001;
m = size(y,1);
n = size(s_init,1);
s_tilde = s_init;
s_hat = s_tilde;
nreqin = 6;

%Load auxilliary variables
num_additional = nargin - nreqin;
if num_additional > 0
    if mod(num_additional,2) == 0
        for i = 1:2:num_additional
            eval([varargin{i}, ' = varargin{', num2str(i+1), '};']);
        end
    else
        error('Wrong number of additional arguments');
    end
end

C_dd_inv_pf = inv(C_dd)^.5;
obj_func = @(s,beta) tikh_obj_eval(y,s,L,C_dd_inv_pf,h_function);

objfunc_value = obj_func(s_tilde);
objfunc_new = objfunc_value;

objfunc_pct = (objfunc_value - objfunc_new)/objfunc_value
s_pct = (sum(abs(s_tilde - s_hat)./s_tilde))/size(s_tilde,1)

num_gradevals = 0;
while ((objfunc_pct > tol_objfunc) && (s_pct > tol_s) && (num_gradevals < max_gradevals)) || (num_gradevals == 0)
    
    format short
    elapsed = etime(clock, start_time);

    s_tilde = s_hat;

    objfunc_value = objfunc_new;
    disp('Iteration         Time Elapsed        Obj. Func');
    disp([num2str(num_gradevals), '      ', num2str(elapsed), '       ', num2str(objfunc_value)]);

    H_tilde = H_tilde_func(s_tilde);
    num_gradevals = num_gradevals + 1;

    h_of_s_tilde = h_function(s_tilde);
    [s_hat] = lin_tikh_inv((y-h_of_s_tilde+H_tilde*s_tilde),C_dd,L,H_tilde);

    objfunc_new = obj_func(s_hat);

    if max_linsearch > 0

        clear linsearch_func
        linsearch_func = @(x) obj_func(s_tilde + (s_hat - s_tilde).*x);
        if objfunc_new < objfunc_value
            disp('Performing linesearch, starting at s_hat');
            lin_options = optimset('MaxIter', max_linsearch,'Display','iter','TolFun',tol_linsearch*objfunc_new);
            step = fminsearch(linsearch_func,1,lin_options)
        else
            disp('Performing linesearch, starting at s_tilde');
            lin_options = optimset('MaxIter', max_linsearch,'Display','iter','TolFun',tol_linsearch*objfunc_value);
            step = fminsearch(linsearch_func,0,lin_options)
        end
        s_hat = s_tilde + step.*(s_hat - s_tilde);
        objfunc_new = obj_func(s_hat,beta_hat);
    end
        
    objfunc_pct = (objfunc_value - objfunc_new)/objfunc_value
    s_pct = max(abs((s_tilde - s_hat)./s_tilde))
    save('qltikh_current_iter','s_hat','H_tilde')
    
end

if objfunc_new < objfunc_value
    s_best = s_hat;
    H_local = H_tilde_func(s_hat);
    objfunc_value = objfunc_new;    
else
    s_best = s_tilde;
    H_local = H_tilde;
end

finish_time = clock;
elapsed = etime(finish_time,start_time);
disp('Iteration         Time Elapsed        Obj. Func');
disp([num2str(num_gradevals), '      ', num2str(elapsed), '       ', num2str(objfunc_value)]);

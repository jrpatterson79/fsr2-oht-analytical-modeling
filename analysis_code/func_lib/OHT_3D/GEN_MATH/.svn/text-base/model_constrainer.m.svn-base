function obj_func = model_constrainer(func_calculator,trial_vec,outside_value,A,b,Aeq,beq,lb,ub,nonlcon)

%model_constrainer: Function which runs (and returns values from) a given
%function ONLY if pre-specified conditions are met. If constraints are not
%met, the function is not run, and the value defined by "outside_value" is
%given.

%[obj_func] =
%model_constrainer(func_calculator,trial_vec,outside_value,A,b,Aeq,beq,lb,ub,nonlcon)
%
%where:
%   -func_calculator is the function that calculates obj_func (as long as
%   constraints are met)
%   -trial_vec is the vector being tried
%   -outside_value is the value that is returend if constraints are not met
%   -A,b,Aeq,beq,lb,ub,and nonlcon are constraints placed on the parameter.
%   See fmincon for more information on supplying these values. These
%   parameters can be passed as empty if not used.

vec_allowed = true;

%Linear / bounds testing cases are simple, all are executed during each
%trial (even if the vector has already failed)
if ~isempty(A)
    if ~isempty(b)
        test1 = ((A*trial_vec) <= b);
        if ~test1 
            vec_allowed = false;
        end
    else
        display('If A is supplied, b must be supplied')
    end
end

if ~isempty(Aeq)
    if ~isempty(beq)
        test2 = ((Aeq*trial_vec) == beq);
        if ~test2
            vec_allowed = false;
        end
    else
        display('If Aeq is supplied, beq must be supplied')
    end
end

if ~isempty(lb)
    test3 = prod(1.*(trial_vec >= lb));
    if ~test3
        vec_allowed = false;
    end
end

if ~isempty(ub)
    test4 = prod(1.*(trial_vec <= ub));
    if ~test4
        vec_allowed = false;
    end
end

%Nonlinear constraint might be time consuming. only tested if other
%constraints have been met.
if vec_allowed == true
    if ~isempty(nonlcon)
        [c, ceq] = nonlcon(trial_vec);
        test5 = prod(1.*(c <= 0));
        test6 = prod(1.*(ceq == 0));
        if (~test5) || (~test6)
            vec_allowed = false;
        end
    end
end

if vec_allowed == true
    obj_func = func_calculator(trial_vec);
else
    obj_func = outside_value;
end
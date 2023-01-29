function [prime_list] = findprime(num_primes)

%findprime: Function that finds a vector of all prime numbers, up to a
%specified length. 
%
%Syntax:
%[prime_list] = findprime(num_primes)
% If num_primes is less than or equal to 1, num_primes returns the first
% prime number, 2.

prime_list = zeros(num_primes,1);

prime_list(1) = 2;
num_found = 1;

trial = 3;
while num_found < num_primes
    %Assume that the trial is prime to begin with
    prime_check = true;
    %Try dividing the trial by each prime found up until now
    for i = 1:1:num_found
        check_val = mod(trial,prime_list(i));
        %If it's perfectly divisible, it's not prime
        if check_val == 0
            prime_check = false;
        end
    end
    if prime_check == true
        prime_list(num_found+1) = trial;
        num_found = num_found+1;
    end
    trial = trial + 1;
end
            
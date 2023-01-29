function [prime_list] = findprime2(max_val)

prime_list = [0 2:1:max_val];

for d = 2:1:ceil(sqrt(max_val))
    prime_list((2*d):d:max_val) = 0;
end
prime_list = prime_list(prime_list>0);




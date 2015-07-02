function [out] = odd(number)

    out = mod(number,2) == 1;

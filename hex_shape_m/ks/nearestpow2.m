function y = nearestpow2(x)
%% najblizsza x potega 2 (rowna lub wyzsza)
a = nextpow2(x); 
%% spr czy x lezy pomiedzy 2^(a-1) i 2^a
if x-2^(a-1) == 2^a-x 
   y = [a-1 a];
%% jesli nie, znajdz najblizsza potege z tych dwoch poteg
else 
   [j,k] = intersect([x-2^(a-1) 2^a-x], min(x-2^(a-1), 2^a-x));
   if k == 1 
      y = a-1;
   else
      y = a;
   end
end
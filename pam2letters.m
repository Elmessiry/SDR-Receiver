% f = pam2letters(seq)
% reconstruct string of +/-1 +/-3 into letters
function f = pam2letters(seq)

S = length(seq);
off = mod(S,4);

if off ~= 0
  sprintf('dropping last %i PAM symbols',off)
  seq = seq(1:S-off);
end

N = length(seq)/4;
for k = 0:N-1
  f(k+1) = base2dec(native2unicode((seq(4*k+1:4*k+4)+99)/2),4);
end

f = native2unicode(f);


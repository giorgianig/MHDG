function res = prolExp(base,exp)
% prolongate the fractional exponential to avoid
% complex number and help NR
% res =  base^exp          -- if base>0
% res = -abs(base)^exp  -- if base<0

res = sign(base).*abs(base).^exp;
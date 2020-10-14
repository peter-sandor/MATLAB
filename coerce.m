function out = coerce(in,lim)
% This function coerces the values of elements of 'in', such that if these are
% beyond certain limits ('in'<minlim or 'in'>maxlim), they will be replaced by
% the appropriate limit.

out=in;
if ~isempty(lim)
    out(in>lim(2))=max(lim);
    out(in<lim(1))=min(lim);
end
end
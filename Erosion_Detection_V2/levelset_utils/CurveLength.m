function l = CurveLength(slice)
[gx, gy] = gradient(slice);
nm = sqrt(gx.^2 + gy.^2);
%nm = nm .* slice;
l = sum(nm(:));
end
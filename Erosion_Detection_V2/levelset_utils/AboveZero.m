function pro_vol = AboveZero(vol)
v_min = min(min(min(vol)));
pro_vol = vol - v_min + 1;
end

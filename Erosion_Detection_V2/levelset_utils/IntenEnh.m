function en_vol = IntenEnh(vol, factor)
vol = single(vol);
vol = AboveZero(vol);
v_max = max(vol(:));
v_min = min(vol(:));
[vx, vy, vz] = gradient(vol);
scale = log(abs(vx) + abs(vy) + abs(vz) + 1) + factor * log(vol + 1);
en_vol = scale .* vol;
en_vol = (en_vol - min(en_vol(:))) / (max(en_vol(:)) - min(en_vol(:)));
en_vol = en_vol * (v_max - v_min) + v_min;
end
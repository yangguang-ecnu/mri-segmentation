function out_vol = add_Rician_vol(vol, nois)

[r, c, s] = size(vol);
out_vol = zeros(r,c,s);
for i=1:s
    out_vol(:,:,i) = add_rician_noise(vol(:,:,i), nois);
end
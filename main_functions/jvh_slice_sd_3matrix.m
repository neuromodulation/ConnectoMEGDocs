function SD=jvh_slice_sd_3matrix(sd_matrix,rois, freq)
arguments
    sd_matrix double
    rois   {mustBeVector}
    freq   {mustBeVector}
end

for r=1:size(sd_matrix,2)
    SD(r,:) =jvh_compute_pooled_sd_multiple_groups_equal_size(squeeze(sd_matrix(r,:,freq)));
end

SD =jvh_compute_pooled_sd_multiple_groups_equal_size(SD(:,rois));

end
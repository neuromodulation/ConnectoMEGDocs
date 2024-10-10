function sd_pooled = jvh_compute_pooled_sd_multiple_groups_equal_size(sd_matrix)
% compute pooled standard deviation of groups with equal variance and size, iterate over 1st dimension,
% sum the sd over the 2nd dimension
arguments
    sd_matrix double
end
repetition_sum = size(sd_matrix,1);
k = size(sd_matrix,2);
sd_pooled = nan(repetition_sum,1);
for i = 1:repetition_sum
sd_pooled(i) = sqrt(sum(sd_matrix(i,:).^2)/k);
end
end

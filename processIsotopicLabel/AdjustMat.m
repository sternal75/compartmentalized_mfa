function r = AdjustMat(d, prob_labeled_carbon_is_C12, prob_natural_C13)
    r = d;
    for i=1:size(d,1)
        v = d(i,:);
        u = adjust(v', length(v)-1, prob_natural_C13, prob_labeled_carbon_is_C12)  ;
        r(i,:) = u';
    end

end
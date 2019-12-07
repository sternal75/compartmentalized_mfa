%file_name = 'input';

[met_list_norm_neg met_name_arr_neg sample_list_short] = ProcessMavenIsotopicLabel(sprintf('data/%s.xls', file_name), 3);
[met_list_norm_pos met_name_arr_pos sample_list_short] = ProcessMavenIsotopicLabel(sprintf('data/%s.xls', file_name), 3);    


if 1==0
    met_list_norm = met_list_norm_neg;
else
    % Merge metabolite lists
    neg_only = setdiff(met_name_arr_neg, met_name_arr_pos);
    pos_only = setdiff(met_name_arr_pos, met_name_arr_neg);
    neg_and_pos = intersect(met_name_arr_pos, met_name_arr_neg);

    for i=1:length(neg_and_pos)
        neg_index = strmatch(neg_and_pos{i}, met_name_arr_neg, 'exact');    
        pos_index = strmatch(neg_and_pos{i}, met_name_arr_pos, 'exact');

        if (length(neg_index)~=1) || (length(pos_index)~=1)
            fprintf('Error..\n');
        end
        if met_list_norm_neg{neg_index}.median_intensity > met_list_norm_pos{pos_index}.median_intensity
            neg_only{end+1} = neg_and_pos{i};
        else
            pos_only{end+1} = neg_and_pos{i};
        end
    end

    met_list_norm = cell(0);
    for i=1:length(neg_only)
        neg_index = strmatch(neg_only{i}, met_name_arr_neg, 'exact');    
        met_list_norm{end+1} = met_list_norm_neg{neg_index};
%         met_list_norm{end}.met_name = sprintf('%s (neg)', met_list_norm{end}.met_name);
    end

    for i=1:length(pos_only)
        pos_index = strmatch(pos_only{i}, met_name_arr_pos, 'exact');    
        met_list_norm{end+1} = met_list_norm_pos{pos_index};
%         met_list_norm{end}.met_name = sprintf('%s (pos)', met_list_norm{end}.met_name);
    end
end





% Fix mas isotopomer distribution vector, by removing natural abundance
% and impurity effects
for i=1:length(met_list_norm)
    v=met_list_norm{i}.data;
    prob_labeled_carbon_is_C12 = 0.01;   % the inpurity of C13
    prob_natural_C13 = 0.011;  %C13 natual abundance
    
    u = AdjustMat(v', prob_labeled_carbon_is_C12, prob_natural_C13);    
    met_list_norm{i}.data  = u';
end





% Find major mass isotopomers
for i=1:length(met_list_norm)
    v = max(met_list_norm{i}.data, [], 2);
    t = find(v >= 0.01);
    met_list_norm{i}.mass_isotopomer_list = t;
end


legend_arr = cell(0);
for x=1:50
    legend_arr{x} = sprintf('m+%d', x-1);
end

mkdir(file_name);
figure;
for i=1:length(met_list_norm)
%    hbar = bar(met_list_norm{i}.data', 'stacked');
    mat_bar = met_list_norm{i}.data(met_list_norm{i}.mass_isotopomer_list, :)';
%     mat_bar(2,:)=0;
    
    mat_bar(end+1,:)=0;
    mat_bar_total =  mat_bar;
   % mat_bar_total = mat_bar .* repmat(met_list_norm{i}.data_total',1,size(mat_bar,2));

    i
%     figure;
    subplot(1,3,i);
    set(gca, 'FontSize', 18);
    bar_label( mat_bar_total, mat_bar);

    
    met_name_for_title = strrep(met_list_norm{i}.met_name,'_',' ');
%     title(met_name_for_title, 'FontSize', 12);
    legend(legend_arr(met_list_norm{i}.mass_isotopomer_list), 'Location', 'NorthEastOutside');

    set(gca, 'XTick', [1:1:length(sample_list_short)]); 
    set(gca, 'XTickLabel', met_name_for_title);
    set(gca, 'Ylim', [0 1]);
    set(gca, 'Xlim', [0.5 1.5]);
%     set(gca, 'FontSize', 6);
    
    s = sprintf('./%s/%s.jpg', file_name, met_list_norm{i}.met_name);
    saveas(gcf, s);
end


function OutputSimulated(idv_opt, EMU, model)
    for(i=1:length(model.mets))
        EMU_indices = find(EMU.list(:,1)==i);
        if(isempty(EMU_indices))
            continue;
        end
        idv_metabolite = idv_opt{EMU_indices(1)};
        idv_metabolite(end+1,:) = 0;
        figure;
        bar_label(idv_metabolite, idv_metabolite);        
        met_name_for_title = strrep(model.mets{i},'_',' ');
        title(met_name_for_title, 'FontSize', 12);
        legend_arr = cell(0);
        for x=1:50
            legend_arr{x} = sprintf('m+%d', x-1);
        end 
        legend(legend_arr(1:size(idv_metabolite,2)), 'Location', 'NorthEastOutside');
%         legend(legend_arr(met_list_norm{i}.mass_isotopomer_list), 'Location', 'NorthEastOutside');

%         set(gca, 'XTick', [1:1:length(sample_list_short)]); 
%         set(gca, 'XTickLabel', sample_list_short);
        set(gca, 'Ylim', [0 1]);
        set(gca, 'Xlim', [0.5 1.5]);
        set(gca, 'FontSize', 6);        
       end

end





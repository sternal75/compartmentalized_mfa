function OutputSimulatedVsMeasured(idv_known, EMU, idv, final_predicted_flux, model, title_output)
    idv_known_arr = zeros(0,1);
    for i=1:length(idv_known) 
        if ~isempty(idv_known{i})
            idv_known_arr = [idv_known_arr; i];
        end
    end
    [idv_opt idv_d cycle_error] = ComputeEmuIDV(EMU, idv, idv_known_arr, final_predicted_flux);
    idv_known_vs_opt=[];
    xLabelForBar=cell(0);
    mass_isotopomer_legend=cell(0);
    for i=1:length(idv_known_arr)
        x = idv_known_arr(i);
        EMU_indices = find(EMU.list(:,1)==x);
        for(j=1:max(cellfun('length',idv_known)))
            mass_isotopomer_legend{end+1}=sprintf('m+%s',num2str(j-1));
        end
        idv_known_vs_opt=[idv_known_vs_opt;zeros(2,max(cellfun('length',idv_known)))];
    %     idv_known_vs_opt=[idv_known_vs_opt;[idv_known{x};idv_opt{EMU_indices(1)}]];
        idv_known_vs_opt(end-1,1:length(idv_known{x}))=idv_known{x};
        idv_known_vs_opt(end,1:length(idv_opt{EMU_indices(1)}))=idv_opt{EMU_indices(1)};
        xLabelForBar{end+1}=strcat(model.mets{x},'_m');
        xLabelForBar{end+1}=strcat(model.mets{x},'_s');
    end
    bar_label(idv_known_vs_opt, idv_known_vs_opt);
    set(gca, 'Xlim', [0 length(idv_known_arr)*2+1]);
    set(gca, 'XTick', [1:1:length(idv_known_arr)*2]);
    set(gca, 'XTickLabel', xLabelForBar);
    legend(mass_isotopomer_legend, 'Location', 'NorthEastOutside');
    %set(gca, 'FontSize', 6);
    set(gca, 'Ylim', [0 1.05]);
    title(title_output);
    drawnow;
end





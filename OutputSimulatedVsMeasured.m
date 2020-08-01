function OutputSimulatedVsMeasured(idv_known, EMU, EMU_met_known_mat, idv, final_predicted_flux, final_predicted_cy_mt_ratio, model, WC_known_metabolites)

    idv_known_vs_opt=zeros(2*length(EMU_met_known_mat),max(cellfun('length',idv)));
    mass_isotopomer_legend=cell(0);        
    for(j=1:max(cellfun('length',idv)))
        mass_isotopomer_legend{end+1}=sprintf('m+%s',num2str(j-1));
    end

    xLabelForBar=cell(0);
    counter=1;
    for i=1:size(EMU_met_known_mat,1)    
        idv_known_vs_opt(counter,1:length(WC_known_metabolites{i}.idv))=WC_known_metabolites{i}.idv;
        counter=counter+1;
        
        x_cy = EMU_met_known_mat(i,1);
        x_mt = EMU_met_known_mat(i,2);

        EMU_indices_cy = find(EMU.list(:,1)==x_cy);
        EMU_indices_mt = find(EMU.list(:,1)==x_mt);
        if(isempty(EMU_indices_cy))
            idv_cy = zeros(size(idv{EMU_indices_mt(1)}));
            iterator_cy_mt_ratio = 0;
        else
            idv_cy = idv{EMU_indices_cy(1)};
        end
        if(isempty(EMU_indices_mt))
            idv_mt = zeros(size(idv{EMU_indices_cy(1)}));
            iterator_cy_mt_ratio = 1;
        else
            idv_mt = idv{EMU_indices_mt(1)};
        end    

        if(((~isnan(x_cy)) && (~isnan(x_mt))))
           iterator_cy_mt_ratio = final_predicted_cy_mt_ratio(i);
        end    
        
        idv_known_vs_opt(counter,1:length(WC_known_metabolites{i}.idv)) = ((iterator_cy_mt_ratio*idv_cy(1:end))+((1-iterator_cy_mt_ratio)*idv_mt(1:end)));
        counter=counter+1;
        counter=counter+1;
        %         idv_known_vs_opt=[idv_known_vs_opt;zeros(2,max(cellfun('length',idv_known)))];
        
        xLabelForBar{end+1}=strcat(strrep(WC_known_metabolites{i}.met_name,'_WC',''),' (m)');
        xLabelForBar{end+1}=strcat(strrep(WC_known_metabolites{i}.met_name,'_WC',''),' (s)');
        xLabelForBar{end+1}='';
    end
    figure;
    bar_label(idv_known_vs_opt, idv_known_vs_opt);
    set(gca, 'Xlim', [0 size(EMU_met_known_mat,1)-1+size(EMU_met_known_mat,1)*2+1]);
    set(gca, 'XTick', [1:1:size(EMU_met_known_mat,1)-1+size(EMU_met_known_mat,1)*2]);
    set(gca, 'XTickLabel', xLabelForBar);
    title('Simulated vs. measured labeling');
    ylabel('Labeling fraction');        
    set(gca, 'FontSize', 22);    
    legend(mass_isotopomer_legend, 'Location', 'NorthEastOutside');
    set(gca, 'Ylim', [0 1.05]);
    set(gcf,'color','w');
    drawnow;    

    
    % plot uptake secretion, measured vs sumulated
    figure;
    xLabelForBar{4}='Serine (m)';
    xLabelForBar{5}='Serine (s)';
    xLabelForBar{6}='';
    xLabelForBar{1}='Glycine (m)';
    xLabelForBar{2}='Glycine (s)';
    xLabelForBar{3}='';
    xLabelForBar{7}='Formate (m)';
    xLabelForBar{8}='Formate (s)';
    
    model.measured_net_fluxes_matrix
    uptake_secretion_known_vs_opt(4) = model.measured_net_fluxes_matrix(1,3);
    uptake_secretion_known_vs_opt(5) = final_predicted_flux(model.measured_net_fluxes_matrix(1,1))-final_predicted_flux(model.measured_net_fluxes_matrix(1,2));
    uptake_secretion_known_vs_opt(6) = 0;
    uptake_secretion_known_vs_opt(1) = -model.measured_net_fluxes_matrix(2,3);
    uptake_secretion_known_vs_opt(2) = final_predicted_flux(model.measured_net_fluxes_matrix(2,2))-final_predicted_flux(model.measured_net_fluxes_matrix(2,1));    
    uptake_secretion_known_vs_opt(3) = 0;
    uptake_secretion_known_vs_opt(7) = -model.measured_net_fluxes_matrix(3,3);
    uptake_secretion_known_vs_opt(8) = -final_predicted_flux(model.measured_net_fluxes_matrix(3,1));        
    bar(uptake_secretion_known_vs_opt);
    set(gca, 'Xlim', [0 9]);
    set(gca, 'XTick', [1:1:8]);
    set(gca, 'XTickLabel', xLabelForBar);
    title('Simulated vs. measured uptake');
    ylabel('Uptake [mM/h]');        
    set(gca, 'FontSize', 22);    
    set(gcf,'color','w');
    drawnow;    
    
end


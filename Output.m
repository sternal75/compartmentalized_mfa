function Output(final_predicted_flux, final_predicted_cy_mt_ratio, min_error, model, the_title)
    %OUTPUT Summary of this function goes here
    %   Detailed explanation goes here
    fluxModelImg = imread('modelpic.jpg');
    text_str    = cell(length(final_predicted_flux),1);
    [dR,tR] = xlsread('model_pic_parameters.xlsx');
    position = [];
    arrow_params = [];
    index_in_model_pix_parameters_file = 1;
    for ind=1:length(final_predicted_flux)
        arrow_params=[arrow_params;dR(index_in_model_pix_parameters_file,4:7)];
        text_str{ind} = sprintf('%s=%s', model.rxns{ind}, num2str(final_predicted_flux(ind),'%0.2f'));
        if(~isempty(findstr(model.rxns{ind},'f')))
            position=[position;dR(index_in_model_pix_parameters_file,1:2)];
        elseif(~isempty(findstr(model.rxns{ind},'b')))
            position=[position;dR(index_in_model_pix_parameters_file,1) dR(index_in_model_pix_parameters_file,2)+11];
            index_in_model_pix_parameters_file = index_in_model_pix_parameters_file+1;
        else
            position=[position;dR(index_in_model_pix_parameters_file,1:2)];
            index_in_model_pix_parameters_file = index_in_model_pix_parameters_file+1;
        end
    end
    
    RGB = insertText(fluxModelImg,position,text_str,'FontSize',12,'BoxColor','white','BoxOpacity',0,'TextColor','black');
    figure;
    opengl hardware;
    %rmse_str{1} = sprintf('RMSE=%s',num2str(sqrt(min_error),'%0.2f'));
%     rmse_str{1} = sprintf('Score=%s',num2str((min_error),'%0.4f'));
%     position = [0 0];
%     RGB = insertText(RGB,position,rmse_str,'FontSize',10,'BoxColor','green','BoxOpacity',0.4,'TextColor','black');
%     position = [0 20];
%     known_met_WC_index=0;
%     for i=1:length(WC_known_metabolites)
%         position(1,2)=position(1,2)+20*known_met_WC_index;
%         x_cy = EMU_met_known_mat(i,1);
%         x_mt = EMU_met_known_mat(i,2);
%         if(((~isnan(x_cy)) && (~isnan(x_mt))))
%            known_met_WC_index = known_met_WC_index+1;
%            cy_mt_ratio_str{1} = sprintf('Cytosol ratio: %s = %.7f\n', WC_known_metabolites{i}.met_name(1:end-3), final_predicted_cy_mt_ratio(i));
%            RGB = insertText(RGB,position,cy_mt_ratio_str,'FontSize',10,'BoxColor','green','BoxOpacity',0.4,'TextColor','black');
%         end        
%     end        
    imshow(RGB);
    hold on
    quiver(arrow_params(:,1),arrow_params(:,2),arrow_params(:,3),arrow_params(:,4),'color','blue', 'MaxHeadSize', 1,'AutoScale','off', 'LineWidth',2.2);
%     quiver(700,500,100,100,'color','red');
    hold off    
    clear title xlabel ylabel;
    title(the_title);

end


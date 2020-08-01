function output_sensitivity_analysis(final_predicted_flux, min_error, model, the_title)
    %OUTPUT Summary of this function goes here
    %   Detailed explanation goes here
    fluxModelImg = imread('modelpic.jpg');
    text_str=cell(0);
    [dR,tR] = xlsread('model_pic_parameters.xlsx');    
    index_in_model_pix_parameters_file = 1;
    index_in_model_rxns = 1;
    position = []; 
    for(i=1:length(final_predicted_flux))
        if(~isnan(final_predicted_flux(i,1)))
            %text_str{end+1} = sprintf('%s=%s-%s',model.rxns{i},num2str(final_predicted_flux(i,1),'%0.0f'),num2str(final_predicted_flux(i,3),'%0.0f'));
            text_str{end+1} = sprintf('%s=[%s %s]',model.rxns{i},num2str(final_predicted_flux(i,1),'%0.0f'),num2str(final_predicted_flux(i,2),'%0.0f'));
        else
%             text_str{end+1} = sprintf('%s=%s',model.rxns{i},num2str(final_predicted_flux(i,2),'%0.0f'));
%             text_str{end+1} = '';
            continue;
        end
        position=[position;dR(index_in_model_pix_parameters_file,:)];
        index_in_model_pix_parameters_file = index_in_model_pix_parameters_file+1;

    end 
%     text_str = cell(length(final_predicted_flux),1);
    for ind=1:length(final_predicted_flux)
%         text_str{ind} = sprintf('%s=%s', model.rxns{ind}, num2str(final_predicted_flux(ind),'%0.2f'));
    end
%     position = [1015 224;1215 266;1215 276;986 466;795 545;795 555;757 515;757 525;838 373;838 383;680 218;680 228;877 198;840 91;840 101;877 147;690 313;690 323;542 313;542 323;393 313;393 323;607 370;505 365;470 140;470 150;479 88;479 98;489 32;489 42;520 390;520 400;314 362;314 372;297 505;297 515;225 473;225 483;128 198;128 208;99 284;325 128;325 138;250 102;149 140;240 32;3 250;3 260;145 597;431 602;500 535;500 545;1138 546;75 416;28 532;145 677;457 464];
    RGB = insertText(fluxModelImg,position,text_str,'FontSize',9,'BoxColor','white','BoxOpacity',0,'TextColor','black');
%     figure;
    opengl hardware;
    %rmse_str{1} = sprintf('RMSE=%s',num2str(sqrt(min_error),'%0.2f'));
    rmse_str{1} = sprintf('Score=%s',num2str((min_error),'%0.4f'));
    position = [5 5];
    RGB = insertText(RGB,position,rmse_str,'FontSize',9,'BoxColor','green','BoxOpacity',0.4,'TextColor','black');
    imshow(RGB);
    clear title xlabel ylabel;
    title(the_title);

end


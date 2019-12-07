%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A computational method for inferring cytosolic and mitochondrial 
% specific metabolic fluxes based on whole-cell level measurements of 
% metabolite isotopic labeling.
% Extending upon standard MFA methodologies, our method further accounts 
% for simulated compartment specific labeling patterns by de-convoluting 
% whole-cell level metabolite isotopic labeling measurements into cytosolic
% and mitochondrial-specific simulated labeling patterns
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result_ratio_best_fit result_ratio_confidece_intervals] = compute_SHMT1_vs_SHMT2_confidence_intervals(experiment_type)
if(experiment_type=='LF')
    file_name = 'input_low_FA';
else
    file_name = 'input_high_FA';
end

addpath(genpath('./cobratoolbox'))

% load metabolic model
LoadModel

% load measured isotope labeling
run processIsotopicLabel/calc;

 
[n,m] = size(model.S);

% convert atom mapping format of test model to a one used by FindEMU
m = [];
for i=1:length(model.rxns)  
    model.mappings_carbon{i}.mapping_r(model.mappings_carbon{i}.mapping_p) = [1:length(model.mappings_carbon{i}.mapping_p)]';
    
    model.mappings_carbon{i}.mapping_mat_p = CreateAtomMappingMat(model, model.mappings_carbon{i}, i);

    m.graph_p = model.mappings_carbon{i}.graph_r;
    m.graph_r = model.mappings_carbon{i}.graph_p;  
    m.mapping_p = model.mappings_carbon{i}.mapping_r;  
    m.mapping_r = model.mappings_carbon{i}.mapping_p;  
    model.mappings_carbon{i}.mapping_mat_r = CreateAtomMappingMat(model, m, i);
end



% Change net flux to foward backward flux - add backward flux
%keep original net reactions model
model_net_fluxes = model;
model_net_fluxes.is_net_flux = zeros(model_net_fluxes.rxn_num,1); %mark net fluxes with '1' and one direction fluxes with 0;
added_fluxes = 0;
measured_net_fluxes_matrix=[];
for i=1:length(model.rxns)
    ind = i+added_fluxes;
    if(~isempty(findstr(model.rxns{ind},'f')))
        model_net_fluxes.is_net_flux(i) = 1;
        backward_flux_name = strrep(model.rxns{ind},'f','b');
        model.rxns = {model.rxns{1:ind} backward_flux_name model.rxns{ind+1:end}};
        model.lb = [model.lb(1:ind); model.lb(ind); model.lb(ind+1:end)];
        model.ub = [model.ub(1:ind); model.ub(ind); model.ub(ind+1:end)];
        model.S = [model.S(:,1:ind) -model.S(:,ind) model.S(:,ind+1:end)];
        model.equality_constraints = [model.equality_constraints(:,1:ind) zeros(size(model.equality_constraints(:,ind))) model.equality_constraints(:,ind+1:end)];
        
        model.mappings_carbon = {model.mappings_carbon{1:ind} model.mappings_carbon{ind} model.mappings_carbon{ind+1:end}};
        model.mappings_carbon{ind+1}.graph_r = model.mappings_carbon{ind}.graph_p;
        model.mappings_carbon{ind+1}.graph_p = model.mappings_carbon{ind}.graph_r; 
        model.mappings_carbon{ind+1}.mapping_r = model.mappings_carbon{ind}.mapping_p';
        model.mappings_carbon{ind+1}.mapping_p = model.mappings_carbon{ind}.mapping_r'; 
        model.mappings_carbon{ind+1}.mapping_mat_r = model.mappings_carbon{ind}.mapping_mat_p;
        model.mappings_carbon{ind+1}.mapping_mat_p = model.mappings_carbon{ind}.mapping_mat_r;
       
        added_fluxes = added_fluxes+1;       
    end

    % build matrix for measured net fluxes - this will be used in the
    % optiomization as part of the objective function
    if(~isempty(model.measured_net_fluxes{i}))
        measured_net_fluxes_matrix(end+1,1)=ind;
        if(~isempty(findstr(model.rxns{ind},'f')))
            measured_net_fluxes_matrix(end,2)=ind+1;
        end
        if(experiment_type=='LF')
            measured_net_fluxes_matrix(end,3)=model.measured_net_fluxes{i}.LF.mean;
            measured_net_fluxes_matrix(end,4)=model.measured_net_fluxes{i}.LF.std;
        else
            measured_net_fluxes_matrix(end,3)=model.measured_net_fluxes{i}.HF.mean;
            measured_net_fluxes_matrix(end,4)=model.measured_net_fluxes{i}.HF.std;                
        end
    end    
end
model.measured_net_fluxes_matrix = measured_net_fluxes_matrix;


% handle equality constraints after adding backward fluxes
for(i=1:size(model.equality_constraints,1))
    ind = find(model.equality_constraints(i,:));
    % backward flux added. An equality constaint should be added also to
    % the backward flux
    if((ind(2)-ind(1))==2)
        add_row = [0 model.equality_constraints(i,1:end-1)];
        model.equality_constraints = [model.equality_constraints;add_row];
    end
end

% model.equality_constraints = [model.equality_constraints;zeros(size(model.equality_constraints(:,1))) model.equality_constraints(:,1:end-1)];
model.rxn_num = length(model.rxns); %update number of reactions
model.used_reactions_status = [1:length(model.rxns)]';  %?? NOT USED??
model.exchange = zeros(model.rxn_num,1);%[0;0;0;0;0;0;0;0];     %?? NOT USED??

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model.mets        = name of metabolites (short names)
% model.metNames    = name of metabolites (full names)
% model.rxns        = name of reactions (short names)
% model.atom_C_num  = number of carbon atoms in each metabolite
% model.met_extra   = binar vector indicating whether metabolites are
%                     external to the model. External metabolites are not
%                     mass-balanced and their isotopic labeling is
%                     assumed to be known.
% model.used_reactions_status = ?? NOT USED??
% model.exchange              = ?? NOT USED??
% model.S                     = stoichiometric matrix
% model.met_num               = number of metabolites (rows in model.S)
% model.rxn_num               = number of reactions (columns in model.S)
% model.mapping_carbon        = atom mapping (cell array whose length is
%                               the number of reactions); see
%                               EMU_implementation_notes.ppt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% EMU_met_known - the EMUs for which experimental mass-isotopomer
% distribution data is exeprimentally available. We also refer to mass-isotopomer
% distribution by IDV (isotopomer distribution vector)
EMU_met_known = cell(model.met_num,1);
EMU_met_known_mat = zeros(0,2);
WC_known_metabolites=cell(0);
for i=1:model.met_num
    EMU_met_known{i} = zeros(model.atom_C_num(i), 0);
end
number_of_measured_wc_idv = 0;
extra_met_isotopomers = cell(model.met_num,1);
for(i=1:length(met_list_norm))
    met_name_iterator = met_list_norm{i}.met_name;
    [num_of_mass_isotopomers num_of_timepoints]=size(met_list_norm{i}.data);
    if(strcmp((met_name_iterator(end-5:end)),'_Media')==1)
        %external metabolite - do nothing
        index=strfind(model.mets,met_name_iterator);
        index=find(~cellfun(@isempty,index));        
        extra_met_isotopomers{index}.atom_idv = met_list_norm{i}.data(:,num_of_timepoints)';
    elseif(strcmp((met_name_iterator(end-2:end)),'_WC')==1)         
        string_rep={'_CY','_MT'};
        EMU_met_known_mat = [EMU_met_known_mat; nan nan];
        WC_known_metabolites{end+1}.idv             = met_list_norm{i}.data(:,num_of_timepoints)'; 
        WC_known_metabolites{end}.idv_variance    = met_list_norm{i}.var(:,num_of_timepoints)';
        WC_known_metabolites{end}.met_name = met_list_norm{i}.met_name;
        for(j=1:length(string_rep))
            met_name_iterator_new = strrep(met_name_iterator,'_WC',string_rep{j});
            index=strfind(model.mets,met_name_iterator_new);
            index=find(~cellfun(@isempty,index));
            if(~(isempty(index)))
                EMU_met_known{index}=ones(model.atom_C_num(index), 1);
                EMU_met_known_mat(end,j) = index;                
            end
        end
        % if the metabolite could be found in both compartments in the
        % model
        if((sum(isnan(EMU_met_known_mat(end,:))))==0)
            number_of_measured_wc_idv = number_of_measured_wc_idv+1;
        end
    else
        fprintf('ERROR - metabolite name must end with _WC or _Media\n');
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find EMUs by going backwards from the starting EMUs ("C" in the example)
fprintf('Find EMU\n');
[EMU_list EMU_met EMU_reactions] = FindEMU (model, EMU_met_known); 

%EMU_list - the list of identified EMUs (20 EMUs in the example network)
%           Each row corresponds to one EMU. The value in the first column
%           is the metabolite in which the EMU is defined. The value in the
%           second column is a running index of the EMU in that metabolite.
%           In the example network there are 6 EMUs for metabolite C (i.e.
%           6 rows in EMU_list whose value in the first column in 3
%           (representing C).
%
% EMU_met - the specific carbons in each EMU. This is a cell array whose
%           length is the number of metabolites in the network. For
%           metabolite i, EMU_met{i} is a matrix whose columns represent
%           the EMUs of metabolite i and the rows represent the carbons.
%           metabolites (a value of 1 represent that the carbon is part of
%           the EMU (see EMU_implementation_notes.ppt).
%
% EMU_reactions - the reactions that produce each EMU. This is a matrix
%                 with 4 columns. Each row i repreesnts a reaction producing
%                 a certain EMU (referred to as an "EMU reaction")
%         EMU_reactions(i,4) - an index of a reaction in the network
%         EMU_reactions(i,3) - an index of an EMU produced by this reaction
%         EMU_reactions(i,1) - an index (row in EMU_list) of an EMU used as
%                              substrate 
%         EMU_reactions(i,2) - an index (row in EMU_list) of a 2nd EMU used as
%                              substrate (in case the combines two EMUs)
%                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Cluster the EMUs, create DAG of clusters, and prepare matrices for efficient computation of IDVs given flux rates
fprintf('Analyze EMU\n');
EMU = AnalyzeEMU(model, EMU_list, EMU_met, EMU_reactions);


% extra_met_isotopomers{1}.atom_idv = [1;1;1];    %Serine_Media
% extra_met_isotopomers{1}.enrichment = 1; 
% extra_met_isotopomers{3}.atom_idv = [0;0;0];    %Glucose
% extra_met_isotopomers{3}.enrichment = 1; 
% extra_met_isotopomers{11}.atom_idv = [0];       %Glycine_Media
% extra_met_isotopomers{11}.enrichment = 1; 
% extra_met_isotopomers{14}.atom_idv = [0];       %Other_1B media
% extra_met_isotopomers{14}.enrichment = 1; 
% % added for the left biderrwectional flux
% extra_met_isotopomers{7}.atom_idv = [0];        %Other_1A media
% extra_met_isotopomers{7}.enrichment = 1; 

extra_met_isotopomers = cell(model.met_num,1);
for(i=1:model.met_num)
    extra_met_isotopomers{i}  = model.met_extra_labeling{i};
end



idv = CreateIDV(EMU, extra_met_isotopomers);

 


min_error=inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run optimization 100 times to find global minima
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_matrix=[];
exitflag_matrix=[];
predicted_flux_matrix=[];
predicted_cy_ratio_matrix=[];
parfor(j=1:100)
    initial_fluxes = rand(length(model.rxns),1).*(model.ub-model.lb)+model.lb;  % Initial flux vector for non-convex optimization
    initial_cy_mt_ratio = rand(number_of_measured_wc_idv,1);
    
    [exitflag error predicted_flux predicted_cy_mt_ratio idv_opt]= ComputeEMUOptFlux(model, EMU, idv, EMU_met_known_mat, met_list_norm, WC_known_metabolites, initial_fluxes, initial_cy_mt_ratio);    
    exitflag_matrix(j) = exitflag;
    if((exitflag ~= -2) && (error < min_error))
        error_matrix(j) = error;
        predicted_flux_matrix(:,j) = predicted_flux;
        predicted_cy_ratio_matrix(:,j) = predicted_cy_mt_ratio;
    end
end
[min_error min_error_index]= min(error_matrix);
best_predicted_flux     = predicted_flux_matrix(:,min_error_index);
best_predicted_cy_ratio = predicted_cy_ratio_matrix(:,min_error_index);
initial_fluxes = best_predicted_flux;
initial_cy_mt_ratio = best_predicted_cy_ratio;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sensitivity analysis on all fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error_array=[]; 
exitflag_array=[];
predicted_fluxes_array=[];
predicted_cy_mt_ratio_array=[];
model.equality_constraints = [model.equality_constraints;zeros(1,size(model.equality_constraints,2))];
[equality_constraints_row equality_constraints_column] = size(model.equality_constraints);

ratio_array = [0.01 0.05 0.1 0.2 0.3 0.4 0.7 1 2 4 7 11 16 22 27 35 45 55 65 75 85 100 120 140 170 200 250 300 400 500 700 1000 1500 2000];
% chi square distribution  - one degree of freedom
CONSTANT_VALUE_FOR_CONFIDENCE_INTERVAL = 3.84;

parfor i=1:length(ratio_array)
    net_flux_ratio = ratio_array(i);
    current_model = model;
    % Sensitivity analysis for SHMT1/SHMT2
    current_model.equality_constraints(end,12)=1;
    current_model.equality_constraints(end,13)=-1;
    current_model.equality_constraints(end,6)=-net_flux_ratio;
    current_model.equality_constraints(end,7)=net_flux_ratio;
%     current_model.equality_constraints(end,27)=net_flux_ratio;
    [exitflag error predicted_flux predicted_cy_mt_ratio idv_opt]= ComputeEMUOptFlux(current_model, EMU, idv, EMU_met_known_mat, met_list_norm, WC_known_metabolites, initial_fluxes, initial_cy_mt_ratio);

    fprintf('\n\t*****************  %f   \t****************\n', net_flux_ratio);    

    error_array(i)=error;
    exitflag_array(i)=exitflag;
    predicted_fluxes_array(:,i)=predicted_flux;
    predicted_cy_mt_ratio_array(:,i)=predicted_cy_mt_ratio;
end
LF_left_right_BD.predicted_fluxes_array=predicted_fluxes_array;
LF_left_right_BD.predicted_cy_mt_ratio_array=predicted_cy_mt_ratio_array;
LF_left_right_BD.error_array=error_array;
LF_left_right_BD.exitflag_array=exitflag_array;
LF_left_right_BD.minimum_error=min_error;

figure;
x=ratio_array;
error_array_within_confidence_intervals = error_array;
error_array_within_confidence_intervals(error_array_within_confidence_intervals>(min(error_array_within_confidence_intervals)+CONSTANT_VALUE_FOR_CONFIDENCE_INTERVAL))=inf;
plot(x(exitflag_array~=-2),error_array(exitflag_array~=-2));
set(gca, 'Ylim', [0 100]);
xlabel('SHMT1/SHMT2','FontSize',30);
ylabel('Score','FontSize',30);
title('Sensitivity analysis SHMT1/SHMT2','FontSize',38);
set(gca, 'FontSize', 30);
set(gca,'xtick',[0:0.5:40])
set(gca, 'Xlim', [0 35]);

result_ratios_within_confience_intervals = ratio_array(error_array_within_confidence_intervals<(min(error_array)+CONSTANT_VALUE_FOR_CONFIDENCE_INTERVAL));
result_ratio_best_fit = (best_predicted_flux(12)-best_predicted_flux(13))/((best_predicted_flux(6)-best_predicted_flux(7)));
result_ratio_confidece_intervals = [min(result_ratios_within_confience_intervals) max(result_ratios_within_confience_intervals)]






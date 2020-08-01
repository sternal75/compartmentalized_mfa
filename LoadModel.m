% clear all;
close all;

% readCbModel - reads input model file in SBML or excell format 
% cobratoolbox
cobraToolBoxModel=readCbModel('input.xlsx');
[dR,tR] = xlsread('input.xlsx');
[dM,tM] = xlsread('input.xlsx','Metabolite List')


% translate to the EMU code model
model.mets      = strrep(cobraToolBoxModel.mets,'[c]','')';
model.metNames  = model.mets;
model.rxns      = cobraToolBoxModel.rxns';
model.lb      = cobraToolBoxModel.lb';
model.ub      = cobraToolBoxModel.ub';
model.S         = full(cobraToolBoxModel.S);
[model.met_num model.rxn_num] = size(model.S);
model.used_reactions_status = [1:length(model.rxns)]';  %?? NOT USED??
model.exchange = zeros(model.rxn_num,1);%[0;0;0;0;0;0;0;0];     %?? NOT USED??
model.equality_constraints = [];
model.non_equality_constraint_values = [];
model.non_equality_constraint_fluxes = [];
% model.forward_backward_flux = dR(:,11);
% model.forward_backward_flux(isnan(model.forward_backward_flux)) = Inf;


model.atom_C_num=[];
model.met_extra=[];
for(i=1:length(model.mets))
    % fill number of carbons for reaction
    idx=strfind(tR(:,3),model.mets{i});
    idx = ~cellfun('isempty',idx);
    met_line_number = find(idx);
    metabolites_in_reaction = strtrim(strsplit(tR{met_line_number(1),3},{' => ',' + '}));
    met_index_in_reaction = strmatch(model.mets(i),metabolites_in_reaction);
    carbon_mapping_in_reaction = strtrim(strsplit(tR{met_line_number(1),14},{' => ',' + '}));
    number_of_carbons = length(carbon_mapping_in_reaction{met_index_in_reaction(1)});
    model.atom_C_num=[model.atom_C_num;number_of_carbons];
    % fill external metabolite = 1
    idx=strfind(tM(:,1),model.mets{i});
    idx = ~cellfun('isempty',idx);
    met_line_number = find(idx);
    if(strcmp(tM(met_line_number,5),'Extra-organism'))
        model.met_extra=[model.met_extra;1];
        % extra-organism labeling
        external_labeling = strsplit(tM{met_line_number,6},'/');
        model.met_extra_labeling{i}.atom_idv=str2num(external_labeling{1});
        model.met_extra_labeling{i}.enrichment=str2num(external_labeling{2});
    else
        model.met_extra=[model.met_extra;0];
        model.met_extra_labeling{i} = [];
    end   
end

model.mappings_carbon = cell(0);
model.measured_net_fluxes=cell(0);
for(i=1:length(model.rxns))
    % init all matrices for carbon mapping
    model.mappings_carbon{end+1}.graph_r.node_info=[];
    model.mappings_carbon{end}.graph_p.node_info=[];
    model.mappings_carbon{end}.graph_r.mets=[];
    model.mappings_carbon{end}.graph_p.mets=[];
    model.mappings_carbon{end}.mapping_p=[];
    temp_reactants =[];
    temp_products  =[];
    % fill carbon mapping
    reactants_and_products = strsplit(tR{i+1,3},{' => '});
    reactants = strtrim(strsplit(reactants_and_products{1},{' + '}));
    products  = strtrim(strsplit(reactants_and_products{2},{' + '}));
    reactants_and_products_carbon_map = strsplit(tR{i+1,14},{' => '});
    reactants_carbon_map = strtrim(strsplit(reactants_and_products_carbon_map{1},{' + '}));
    products_carbon_map  = strtrim(strsplit(reactants_and_products_carbon_map{2},{' + '}));    
     
    % add to stoichiometric matrix - to handle chiral metabolites. for e.g.
    % fumarate(1234)==>malate(1234) and fumarate(4321)==>malate(1234),
    % which I added to the model as two reactions for simplicity
    addDuplicateReactionsToSMatrix=zeros(1,length(model.rxns));
    for(j=i+1:length(model.rxns))
        %if(strcmp(tR(j+1,3),tR(i+1,3)))
        if((~isnan(dR(i,10))) & (dR(j,10)==dR(i,10))) 
            addDuplicateReactionsToSMatrix(i)=1;
            addDuplicateReactionsToSMatrix(j)=-1;
            model.equality_constraints = [model.equality_constraints;addDuplicateReactionsToSMatrix];
            break;
        end
    end 
    
    % go over all reactants
    for(j=1:length(reactants))
        idx = strfind(model.mets,reactants{j});
        idx = ~cellfun('isempty',idx);
        met_index = find(idx);
        model.mappings_carbon{end}.graph_r.mets = [model.mappings_carbon{end}.graph_r.mets;met_index];
        model.mappings_carbon{end}.graph_r.node_info(end+1:end+model.atom_C_num(idx),1)=j;
        
        if(j==1)
            reactant_mapping=str2double(regexp(num2str(reactants_carbon_map{j}),'\d','match'));
        else
            reactant_mapping=[reactant_mapping str2double(regexp(num2str(reactants_carbon_map{j}),'\d','match'))+reactant_mapping(end)];
        end        
    end
    % go over all products
    for(j=1:length(products))
        idx = strfind(model.mets,products{j});
        idx = ~cellfun('isempty',idx);
        met_index = find(idx);       
        model.mappings_carbon{end}.graph_p.mets = [model.mappings_carbon{end}.graph_p.mets;met_index];
        model.mappings_carbon{end}.graph_p.node_info(end+1:end+model.atom_C_num(idx),1)=j;
        
        if(j==1)
            product_mapping=str2double(regexp(num2str(products_carbon_map{j}),'\d','match'));
        else
            product_mapping=[product_mapping str2double(regexp(num2str(products_carbon_map{j}),'\d','match'))+product_mapping(end)];
        end        
    end    
    model.mappings_carbon{end}.mapping_p(product_mapping,1)=reactant_mapping';
    
    % load measured fluxes
    
    if(~isnan(dR(i,11)))
        model.measured_net_fluxes{end+1}.HF.mean    = dR(i,11);
        model.measured_net_fluxes{end}.HF.std     = dR(i,12);
        model.measured_net_fluxes{end}.LF.mean    = dR(i,13);
        model.measured_net_fluxes{end}.LF.std     = dR(i,14);
    else
        model.measured_net_fluxes{end+1}=[];
    end
    
end


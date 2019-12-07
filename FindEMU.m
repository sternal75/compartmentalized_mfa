function [EMU_list EMU_met EMU_reactions] = FindEMU (model, EMU_met_param)
global EMU_list EMU_met EMU_reactions EMU_list_ind EMU_hash private_EMU_hash_vec fh counter;

fh = fopen('FindIsotopomers.txt', 'wt');
% fh=1;

% [met_num, rxn_num] = size(model.S);
met_num = length(model.mets);
rxn_num = length(model.rxns);
EMU_list = [];
EMU_list_ind = 1;
EMU_met = EMU_met_param;
EMU_reactions = zeros(0,4);
EMU_hash = cell(met_num,1);
private_EMU_hash_vec = 2.^[0:100];
private_EMU_hash_vec = private_EMU_hash_vec';
counter = 0;

for i=1:met_num
    EMU_hash{i} = zeros(0,2);
    for x=1:size(EMU_met{i},2)
        AddIsotopomerHash(i, EMU_met{i});
    end
end


S = model.S;
S(:, setdiff([1:rxn_num], model.used_reactions_status)) = 0;

while (EMU_list_ind <= size(EMU_list,1))
    if mod(counter, 1000) == 0
        fprintf('%d\n', counter);
    end
    
    target_met = EMU_list(EMU_list_ind, 1);
    target_met_EMU_ind = EMU_list(EMU_list_ind, 2);
    
    if model.met_extra(target_met)
        EMU_list_ind = EMU_list_ind+1;
        continue;
    end
    
    target_EMU = EMU_met{target_met}(:, target_met_EMU_ind);
    s = DispEMU(model.metNames{target_met}, target_EMU); 
    fprintf(fh,'\n%s\n',s);

    p = find(S(target_met, :)>0);
    for x = p
        m = model.mappings_carbon{x};            
        AddIsotopomers(model, target_met, target_met_EMU_ind, target_EMU, m, x, 1);
    end

    EMU_list_ind = EMU_list_ind+1;
end
    
%fclose(fh);





function AddIsotopomers(model, target_met, target_met_EMU_ind, target_EMU, m, rec_num, is_reverse, forward)
global EMU_list EMU_met EMU_reactions EMU_list_ind EMU_hash private_EMU_hash_vec fh counter;

[met_num, rxn_num] = size(model.S);

EMU_reactions(end+1, :) = [0, 0, EMU_list_ind, rec_num];

% Single mapping
t = find(m.mapping_mat_r.single.met == target_met);
EMU_reactions_pos = 1;
source_EMU_counter = 0;
if (~isempty(t))    
    for z = 1:length(t)
        y = t(z);
        source_EMU = m.mapping_mat_r.single.mat{y} * target_EMU;
        if sum(source_EMU) == 0
            continue;
        end
        source_met = m.graph_r.mets(y);
        [found source_EMU_ind] = AddIsotopomerHash(source_met, source_EMU);
        EMU_reactions(end, EMU_reactions_pos) = source_EMU_ind;
        EMU_reactions_pos = EMU_reactions_pos+1;

        if (found == 0)
            EMU_met{source_met}(:, end+1) = source_EMU;
            counter = counter+1;
        end
        
        s = DispEMU(model.metNames{source_met}, source_EMU);    
        if(source_EMU_counter==0)
            fprintf(fh,'\t\t<- %s', s);
        else
            fprintf(fh,' + %s', s);
        end
        source_EMU_counter = source_EMU_counter+1;
    end
    fprintf(fh, '\n');
end

% Double mapping
[t1,t2] = find(m.mapping_mat_r.double.met == target_met);        
for y = [t1, t2]'
    target_EMU_ind = y(1);
    
    source_EMU = m.mapping_mat_r.double.mat{y(2), target_EMU_ind} * target_EMU;
    if sum(source_EMU) == 0
        continue;
    end
    source_met = m.graph_r.mets(y(2));    
    
    [found source_EMU_ind] = AddIsotopomerHash(source_met, source_EMU);
    EMU_reactions(end, 1) = source_EMU_ind;
        
    if (found == 0)
        EMU_met{source_met}(:, end+1) = source_EMU;
        counter = counter+1;
    end        
    s = DispEMU(model.metNames{source_met}, source_EMU);      
    fprintf(fh,'\t\t<- %s', s);
    fprintf(fh, '\n');

end



function [found isotopomer_index] = AddIsotopomerHash(met, isotopomer)
global EMU_hash private_EMU_hash_vec EMU_list;

n = length(isotopomer);

n = sum(private_EMU_hash_vec(1:n).*isotopomer);
t = find(EMU_hash{met}(:,1) == n);
if isempty(t)
    EMU_hash{met}(end+1,:) = [n size(EMU_list,1)+1];
    EMU_list = [EMU_list; met size(EMU_hash{met},1) n];
    isotopomer_index = size(EMU_list,1);
    found = 0;
else
    isotopomer_index = EMU_hash{met}(t,2);
    found = 1;
end



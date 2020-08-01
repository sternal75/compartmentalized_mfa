function EMU = AnalyzeEMU(model, EMU_list, EMU_met, EMU_reactions)
[model.met_num model.rxn_num] = size(model.S);

met_emus = zeros(model.met_num,1);
met_vars = zeros(model.met_num,1);
EMU_name = cell(0);
for i=1:length(EMU_list)
    met = EMU_list(i,1);
    met_isotope = EMU_list(i,2);
    
    met_emus(met) = size(EMU_met{met},2);

    EMU_size(i) = sum(EMU_met{met}(:, met_isotope));
    met_vars(met) = met_vars(met) + EMU_size(i)+1;
    fprintf('%d:\t', i);
    EMU_name{i} = DispEMU(model.mets{met}, EMU_met{met}(:, met_isotope)); 
    fprintf('%s\n', EMU_name{i});
end
fprintf('Mets with EMU = %d (out of %d)\n', sum(met_emus~=0), length(model.mets));
fprintf('Total EMU = %d\n', sum(met_emus));
fprintf('Total vars = %d\n', sum(met_vars));



% Create isotope interactions map
n = length(EMU_list);

EMU_mat_1st_reaction = sparse(EMU_reactions(:, 1), EMU_reactions(:,3), ones(length(EMU_reactions),1), n,n);
EMU_mat_2st_reaction = sparse(EMU_reactions(find(EMU_reactions(:,2)~=0),2), EMU_reactions(find(EMU_reactions(:,2)~=0),3), ones(length(find(EMU_reactions(:,2)~=0)),1), n,n);
EMU_mat=EMU_mat_1st_reaction+EMU_mat_2st_reaction;

% fix size of EMU_mat
m = EMU_mat;
[i,j,s] = find(m);
m = sparse(i,j,s,n,n);
EMU_mat = abs(m);

% find connected components
[s,c]=graphconncomp(m);
%[a,b] = sort(c);
%m = m(b,b);

%v = find (c==57);
%DispEMUs(model, EMU_met, EMU_list, v);
%return;
    
% reduce components
m1 = sparse(s, n);
for i=1:s
    v = sum(m(c==i, :), 1);
    m1(i, :) = v;
end
m2 = sparse(s, s);
for i=1:s
    v = sum(m1(:, c==i), 2);
    m2(:, i) = v;
end
m = m2;
m(m~=0) = 1;

% create DAG
v = diag(m)+0;
m=m-diag(v);
o = graphtopoorder(sparse(m));  
% view(biograph(sparse(m)))

EMU_cluster_order = o;

for i=1:length(o)
    x = o(i);
    v = find (c==x);
    fprintf('Cluster #%d:\t', x);
    DispEMUs(model, EMU_met, EMU_list, v); 
end

% Create EMU matrices

EMU_cluster = cell(s,1);
for i=1:s
    EMU_cluster{i}.EMU_ind = find(c==i);            % EMU indecies in cluster
    EMU_cluster{i}.EMU_num = length(EMU_cluster{i}.EMU_ind); % number of EMUs in cluster
    EMU_cluster{i}.EMU_cluster_ind = sparse(n,1);
    EMU_cluster{i}.EMU_cluster_ind(EMU_cluster{i}.EMU_ind) = [1:EMU_cluster{i}.EMU_num]'; 
    
    %EMU_cluster{i}.prev_EMU = find(sum(abs(EMU_mat(:, EMU_cluster{i}.EMU_ind)),2));
    %EMU_cluster{i}.prev_EMU = setdiff(EMU_cluster{i}.prev_EMU, EMU_cluster{i}.EMU_ind)';
    %EMU_cluster{i}.prev_EMU_num = length(EMU_cluster{i}.prev_EMU);
    
    t1 = zeros(0,1); t2 = zeros(0,1);
    for x=EMU_cluster{i}.EMU_ind
        t1 = [t1; find(EMU_reactions(:,3) == x & EMU_reactions(:,2) == 0)];
        t2 = [t2; find(EMU_reactions(:,3) == x & EMU_reactions(:,2) ~= 0)];
    end
    
    %Reactions with one EMU entering the cluster
    EMU_cluster{i}.prev_EMU_ind = unique(EMU_reactions(t1,1));
    EMU_cluster{i}.prev_EMU_ind = setdiff(EMU_cluster{i}.prev_EMU_ind, EMU_cluster{i}.EMU_ind)';
    EMU_cluster{i}.prev_EMU_num = length(EMU_cluster{i}.prev_EMU_ind);
    EMU_cluster{i}.prev_EMU_cluster_ind = sparse(n,1);
    EMU_cluster{i}.prev_EMU_cluster_ind(EMU_cluster{i}.prev_EMU_ind) = [1:EMU_cluster{i}.prev_EMU_num]';
        
    %Reactions with two EMU entering the cluster
    EMU_cluster{i}.prev_EMU_pair_ind = unique(EMU_reactions(t2,1:2), 'rows');
    EMU_cluster{i}.prev_EMU_pair_num = size(EMU_cluster{i}.prev_EMU_pair_ind,1);
    EMU_cluster{i}.prev_EMU_pair_cluster_ind = sparse(EMU_cluster{i}.prev_EMU_pair_ind(:,1), EMU_cluster{i}.prev_EMU_pair_ind(:,2), [1:EMU_cluster{i}.prev_EMU_pair_num]',n,n);
    
    EMU_cluster{i}.left_mat = zeros(0,3);
    EMU_cluster{i}.right_mat1 = zeros(0,3);
    EMU_cluster{i}.right_mat2 = zeros(0,3);
end

for i=1:size(EMU_reactions,1)
    source_EMU = EMU_reactions(i,1);
    source_EMU_cluster = c(source_EMU);

    target_EMU = EMU_reactions(i,3);
    target_EMU_cluster = c(target_EMU);

    reaction_EMU = EMU_reactions(i,4);
                                                    
        EMU_cluster{target_EMU_cluster}.left_mat = [EMU_cluster{target_EMU_cluster}.left_mat;...
                                                    EMU_cluster{target_EMU_cluster}.EMU_cluster_ind(target_EMU),...
                                                    EMU_cluster{target_EMU_cluster}.EMU_cluster_ind(target_EMU),...
                                                    -reaction_EMU];
                                                
    if (source_EMU_cluster == target_EMU_cluster)
        fprintf('(%d) Cluster #%d: %s -> Cluster #%d %s\n', i, source_EMU_cluster, EMU_name{source_EMU}, target_EMU_cluster, EMU_name{target_EMU});
        EMU_cluster{source_EMU_cluster}.left_mat = [EMU_cluster{source_EMU_cluster}.left_mat;...
                                                    EMU_cluster{source_EMU_cluster}.EMU_cluster_ind(target_EMU),...
                                                    EMU_cluster{source_EMU_cluster}.EMU_cluster_ind(source_EMU),...
                                                    reaction_EMU];

        continue;
    end
    
    if (EMU_reactions(i,2) == 0)
        fprintf('(%d) Cluster #%d: %s -> Cluster #%d %s\n', i, source_EMU_cluster, EMU_name{source_EMU}, target_EMU_cluster, EMU_name{target_EMU});
        EMU_cluster{target_EMU_cluster}.right_mat1 = [EMU_cluster{target_EMU_cluster}.right_mat1;...
                                                    EMU_cluster{target_EMU_cluster}.EMU_cluster_ind(target_EMU),...
                                                    EMU_cluster{target_EMU_cluster}.prev_EMU_cluster_ind(source_EMU),...
                                                    -reaction_EMU];
        continue;
    end
    
    source_EMU2 = EMU_reactions(i,2);
    source_EMU_cluster2 = c(source_EMU2);
    fprintf('(%d) Cluster #%d: %s + Cluster #%d: %s -> Cluster #%d %s\n', i, source_EMU_cluster, EMU_name{source_EMU}, source_EMU_cluster2, EMU_name{source_EMU2}, target_EMU_cluster, EMU_name{target_EMU});
    
    EMU_cluster{target_EMU_cluster}.right_mat2 = [EMU_cluster{target_EMU_cluster}.right_mat2;...
                                                EMU_cluster{target_EMU_cluster}.EMU_cluster_ind(target_EMU),...
                                                EMU_cluster{target_EMU_cluster}.prev_EMU_pair_cluster_ind(source_EMU, source_EMU2),...
                                                -reaction_EMU];
end

%%%%%%%%%%%%%%%%
% DEBUG
fprintf('Left Matrix:\n');
DispMat(model, EMU_name, EMU_cluster{4}.left_mat, EMU_cluster{4}.EMU_num, EMU_cluster{4}.EMU_num);

fprintf('Right Matrix 1:\n');
DispMat(model, EMU_name, EMU_cluster{4}.right_mat1, EMU_cluster{4}.EMU_num, EMU_cluster{4}.prev_EMU_num);

fprintf('Right Matrix 2:\n');
DispMat(model, EMU_name, EMU_cluster{4}.right_mat2, EMU_cluster{4}.EMU_num, EMU_cluster{4}.prev_EMU_pair_num);
%%%%%%%%%%%%%%%%

for i=1:s
    EMU_cluster{i}.left_mat_d = cell(0);
    EMU_cluster{i}.right_mat_d1 = cell(0);
    EMU_cluster{i}.right_mat_d2 = cell(0);
    for x=1:model.rxn_num
        t = find(abs(EMU_cluster{i}.left_mat(:,3)) == x);
        EMU_cluster{i}.left_mat_d{x} = sparse(EMU_cluster{i}.left_mat(t,1), EMU_cluster{i}.left_mat(t,2), sign(EMU_cluster{i}.left_mat(t,3)), EMU_cluster{i}.EMU_num, EMU_cluster{i}.EMU_num);

        t = find(abs(EMU_cluster{i}.right_mat1(:,3)) == x);
        EMU_cluster{i}.right_mat_d1{x} = sparse(EMU_cluster{i}.right_mat1(t,1), EMU_cluster{i}.right_mat1(t,2), sign(EMU_cluster{i}.right_mat1(t,3)), EMU_cluster{i}.EMU_num, EMU_cluster{i}.prev_EMU_num);

        t = find(abs(EMU_cluster{i}.right_mat2(:,3)) == x);
        EMU_cluster{i}.right_mat_d2{x} = sparse(EMU_cluster{i}.right_mat2(t,1), EMU_cluster{i}.right_mat2(t,2), sign(EMU_cluster{i}.right_mat2(t,3)), EMU_cluster{i}.EMU_num, EMU_cluster{i}.prev_EMU_pair_num);
    end
end


EMU.cluster = EMU_cluster;
EMU.list = EMU_list;
EMU.met = EMU_met;
EMU.reactions = EMU_reactions;
EMU.name = EMU_name;
EMU.size = EMU_size;
EMU.cluster_order = EMU_cluster_order;

function DispMat(model, EMU_name, mat, sx, sy)
for x=1:sx
    for y=1:sy
        v = find(mat(:,1) == x & mat(:,2) == y);
        if (length(v) == 0)
            fprintf('0');
        end
        for z=v'
            v = mat(z,3);
            fprintf('%d%s ', sign(v)+0, model.rxns{abs(v)});
        end
        fprintf('\t\t\t');
    end
    fprintf('\n');    
end

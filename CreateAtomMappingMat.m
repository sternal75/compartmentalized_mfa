function mapping_mat_p = CreateAtomMappingMat(model, m, x)
    if (x==434)
        a=11;
    end
    
    mapping_mat_p.single.met = [];
    mapping_mat_p.single.met_ind = [];
    mapping_mat_p.double.met = [];
    mapping_mat_p.double.met_ind = [];
    
    target_met_list = m.graph_p.node_info(:,1);
    source_met_list = m.graph_r.node_info(:,1);
    
    u = [m.graph_r.node_info(m.mapping_p, 1), target_met_list];
    u = unique(u, 'rows');
    for i=u(:,2)'
        t = find(u(:,2) == i);
        if length(t) > 2
            fprintf('Reaction %d has more than 2 mets with carbon mappings:\t',x);
            DispReaction(model, x);
            fprintf('\n');
        end
        
        % Target metabolite mapped to a single source metabolite 
        if length(t) == 1
            source_met_ind = u(t, 1);
            target_met_ind = i;

            source_met = m.graph_r.mets(source_met_ind);
            target_met = m.graph_p.mets(target_met_ind);
            
            source_carbons = model.atom_C_num(source_met);
            target_carbons = model.atom_C_num(target_met);
            
            v_r = find(m.graph_r.node_info(:,1) == source_met_ind);
            v_p = find(m.graph_p.node_info(:,1) == target_met_ind);
            v = m.mapping_p(v_p) - v_r(1) + 1;
            
            mapping_mat_p.single.met_ind(target_met_ind) = source_met_ind;
            mapping_mat_p.single.met(target_met_ind) = source_met;            
            mapping_mat_p.single.mat{target_met_ind} = sparse([1:target_carbons], v, ones(target_carbons,1), target_carbons, source_carbons);
        end
        
        % Target metabolite mapped to two source metabolite 
        if length(t) == 2
            source_met1_ind = u(t(1), 1);
            source_met2_ind = u(t(2), 1);
            target_met_ind = i;

            source_met1 = m.graph_r.mets(source_met1_ind);
            source_met2 = m.graph_r.mets(source_met2_ind);
            target_met = m.graph_p.mets(target_met_ind);
            
            source_carbons1 = model.atom_C_num(source_met1);
            source_carbons2 = model.atom_C_num(source_met2);
            target_carbons = model.atom_C_num(target_met);

            v_r1 = find(m.graph_r.node_info(:,1) == source_met1_ind);
            v_r2 = find(m.graph_r.node_info(:,1) == source_met2_ind);
            v = find(m.graph_p.node_info(:,1) == target_met_ind);
            v_p1 = find(m.graph_p.node_info(:,1) == target_met_ind & m.graph_r.node_info(m.mapping_p,1) == source_met1_ind);
            v_p2 = find(m.graph_p.node_info(:,1) == target_met_ind & m.graph_r.node_info(m.mapping_p,1) == source_met2_ind);
            v1 = zeros(target_carbons, 1);
            v2 = zeros(target_carbons, 1);            
            
%            t = find(m.mapping_p(v_p));
            v1(v_p1-v(1)+1) = m.mapping_p(v_p1) - v_r1(1) + 1;
            
 %           t = find(m.mapping_p(v_p));
            v2(v_p2-v(1)+1) = m.mapping_p(v_p2) - v_r2(1) + 1;
            
            mapping_mat_p.double.met_ind(:, target_met_ind) = [source_met1_ind; source_met1_ind];
            mapping_mat_p.double.met(:, target_met_ind) = [source_met1; source_met2];
            t = find(v1);
            mapping_mat_p.double.mat{target_met_ind,1} = sparse(t, v1(t), ones(length(t),1), target_carbons, source_carbons1);
            t = find(v2);
            mapping_mat_p.double.mat{target_met_ind,2} = sparse(t, v2(t), ones(length(t),1), target_carbons, source_carbons2);
        end
    end

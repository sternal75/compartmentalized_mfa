function [idv result_idv_d cycle_error] = ComputeEmuIDV(EMU, idv, idv_d_emu_ind, flux)

tic_begin = tic;
ticks = 0;

n = length(EMU.list);
cycle_error = 0;
rxn_num = length(flux);

result_idv_d_arr = zeros(n, 1);
for i=1:length(idv_d_emu_ind)
    result_idv_d_arr(find(EMU.list(:,1)==idv_d_emu_ind(i)))=idv_d_emu_ind(i);
end
idv_d_emu_num = length(idv_d_emu_ind);
% result_idv_d_arr = zeros(n, 1);
% result_idv_d_arr(idv_d_emu_ind) = [1:idv_d_emu_num];

% result_idv_d = cell(idv_d_emu_num, 1);
result_idv_d = cell(length(result_idv_d_arr), 1);

for i=idv_d_emu_ind(:)'
%     result_idv_d{i} = sparse(rxn_num, EMU.size(i)+1);
end

idv_d_initial = cell(n, 1);
for i=1:n
    if (~isempty(idv{i}))
        idv_d_initial{i} = zeros(1, length(idv{i}));
    end
end

cluster_num = length(EMU.cluster);
left_mat = cell(cluster_num,1);
right_mat1 = cell(cluster_num,1);
right_mat_idv1 = cell(cluster_num,1);
right_mat2 = cell(cluster_num,1);
right_mat_idv2 = cell(cluster_num,1);

for i=1:cluster_num
   clust_ind = EMU.cluster_order(i);
   
   if (isempty(EMU.cluster{clust_ind}.left_mat))
       %fprintf('Cluster #%d: external EMUs\n', clust_ind);
       continue;
   end
   
   cluster_EMU_size = EMU.size(EMU.cluster{clust_ind}.EMU_ind(1));
   left_mat{i} = sparse(   EMU.cluster{clust_ind}.left_mat(:,1), ...
                        EMU.cluster{clust_ind}.left_mat(:,2), ...
                        flux(abs(EMU.cluster{clust_ind}.left_mat(:,3))) .* sign(EMU.cluster{clust_ind}.left_mat(:,3)), ...
                        EMU.cluster{clust_ind}.EMU_num, EMU.cluster{clust_ind}.EMU_num);
                    
   right_mat1{i} = sparse( EMU.cluster{clust_ind}.right_mat1(:,1), ...
                        EMU.cluster{clust_ind}.right_mat1(:,2), ...
                        flux(abs(EMU.cluster{clust_ind}.right_mat1(:,3))) .* sign(EMU.cluster{clust_ind}.right_mat1(:,3)), ...
                        EMU.cluster{clust_ind}.EMU_num, EMU.cluster{clust_ind}.prev_EMU_num);

    right_mat_idv1{i} = sparse(EMU.cluster{clust_ind}.prev_EMU_num, cluster_EMU_size+1);
    if ~isempty(right_mat1{i})
        for z=1:EMU.cluster{clust_ind}.prev_EMU_num
            right_mat_idv1{i}(z, :) = idv{EMU.cluster{clust_ind}.prev_EMU_ind(z)}';
        end
    end
                    
   right_mat2{i} = sparse( EMU.cluster{clust_ind}.right_mat2(:,1), ...
                        EMU.cluster{clust_ind}.right_mat2(:,2), ...
                        flux(abs(EMU.cluster{clust_ind}.right_mat2(:,3))) .* sign(EMU.cluster{clust_ind}.right_mat2(:,3)), ...
                        EMU.cluster{clust_ind}.EMU_num, EMU.cluster{clust_ind}.prev_EMU_pair_num);
                    
   right_mat_idv2{i} = sparse(EMU.cluster{clust_ind}.prev_EMU_pair_num, cluster_EMU_size+1);
   if ~isempty(right_mat2{i})
       [aa, bb, cc] = find( EMU.cluster{clust_ind}.prev_EMU_pair_cluster_ind );
       for z = 1:length(aa)
           EMU1 = aa(z);
           EMU2 = bb(z);       
           EMU_pair_ind = cc(z);       
           m = CreateCauchyMat( EMU.size(EMU1)+1, EMU.size(EMU2)+1);    % Move to pre-processing
           [a,b,c] = find(m);
           m = sparse(a, b, idv{EMU1}(c));

           v = m*idv{EMU2}';
           right_mat_idv2{i}(EMU_pair_ind, :) = v;
       end
   end
   
   %m1 = inv(left_mat);
   m2 = right_mat1{i}*right_mat_idv1{i};
   m3 = right_mat2{i}*right_mat_idv2{i};
   m4 = m2+m3;
   
   if (1==0)
       max_cluster_in_flux = max(max(abs(m4)));
       max_cluster_flux = max(max(abs(left_mat{i})));   
   
        if (max_cluster_in_flux < 1e-3 & max_cluster_flux > 1e-8)
            %cycle_error = cycle_error + max(log10(max_cluster_flux/max_cluster_in_flux)-3,0);
            cycle_error = 1;
            cycle_error = min(cycle_error, 10);

            for x=1:EMU.cluster{clust_ind}.EMU_num
               idv{EMU.cluster{clust_ind}.EMU_ind(x)} = zeros(1, cluster_EMU_size+1);
               %fprintf('%s = %s\t', EMU.name{EMU.cluster{clust_ind}.EMU_ind(x)}, DispIDV(m(x,:)));
            end
           continue;
        end
   end
   
   %r = rank(left_mat{i}+0); %%% REMOVE ?!
   %a = size(left_mat{i},1);
   %if (r ~= a)
%       z=1;
%   end

   %used_met = find( sum(abs(left_mat),1)>1e-5 );
   t1 = tic;
   m = left_mat{i}\m4;   

   t2 = toc(t1);
   ticks = ticks + t2;
   
   %fprintf('Cluster #%d:\t', clust_ind);
   for x=1:EMU.cluster{clust_ind}.EMU_num
       idv{EMU.cluster{clust_ind}.EMU_ind(x)} = m(x,:);
    if(isnan(idv{EMU.cluster{clust_ind}.EMU_ind(x)}(1)))
        alon=1
    end       
       %fprintf('%s = %s\t', EMU.name{EMU.cluster{clust_ind}.EMU_ind(x)}, DispIDV(m(x,:)));
   end
   %fprintf('\n');
end

%%%%%%%%%%%%%%
% Compute derivatives

for x = 1:rxn_num
    idv_d = idv_d_initial;

    for i=1:length(EMU.cluster)
       clust_ind = EMU.cluster_order(i);
       cluster_EMU_size = EMU.size(EMU.cluster{clust_ind}.EMU_ind(1));  
       if (isempty(EMU.cluster{clust_ind}.left_mat))
           continue;
       end
   
        right_mat_idv_d1 = sparse(EMU.cluster{clust_ind}.prev_EMU_num, cluster_EMU_size+1);
        if ~isempty(right_mat1{i})
            for z=1:EMU.cluster{clust_ind}.prev_EMU_num
                right_mat_idv_d1(z, :) = idv_d{EMU.cluster{clust_ind}.prev_EMU_ind(z)}';
            end
        end

       right_mat_idv_d2 = sparse(EMU.cluster{clust_ind}.prev_EMU_pair_num, cluster_EMU_size+1);
       if ~isempty(right_mat2{i})
           [aa, bb, cc] = find( EMU.cluster{clust_ind}.prev_EMU_pair_cluster_ind );
           for z = 1:length(aa)
               EMU1 = aa(z);
               EMU2 = bb(z);       
               EMU_pair_ind = cc(z);       
               m = CreateCauchyMat( EMU.size(EMU1)+1, EMU.size(EMU2)+1);
               [a,b,c] = find(m);
               
               m = sparse(a, b, idv_d{EMU1}(c));
               v1 = m*idv{EMU2}';

               m = sparse(a, b, idv{EMU1}(c));
               v2 = m*idv_d{EMU2}';
               v = v1+v2;

               right_mat_idv_d2(EMU_pair_ind, :) = v;
           end
       end
       
       left_mat_idv = sparse(EMU.cluster{clust_ind}.EMU_num, cluster_EMU_size+1);
       for z=1:EMU.cluster{clust_ind}.EMU_num
           left_mat_idv(z, :) = idv{EMU.cluster{clust_ind}.EMU_ind(z)};
       end
       
       m1 = EMU.cluster{clust_ind}.left_mat_d{x} * left_mat_idv;
       m2 = EMU.cluster{clust_ind}.right_mat_d1{x}*right_mat_idv1{i} + right_mat1{i}*right_mat_idv_d1;
       m3 = EMU.cluster{clust_ind}.right_mat_d2{x}*right_mat_idv2{i} + right_mat2{i}*right_mat_idv_d2;       
       m4 = m2+m3-m1;
       t1 = tic;
       m = left_mat{i}\m4;       
       t2 = toc(t1);
       ticks = ticks + t2;
       
       for z=1:EMU.cluster{clust_ind}.EMU_num
           emu_ind = EMU.cluster{clust_ind}.EMU_ind(z);
           idv_d{emu_ind} = m(z,:);
           
           if result_idv_d_arr(emu_ind)
%                 result_idv_d{result_idv_d_arr(emu_ind)}(x, :) = m(z,:);
                result_idv_d{emu_ind}(x, :) = m(z,:);
           end
       end

    end
   
end
%ticks
tic_end = toc(tic_begin);
%fprintf('Fraction of time in inv = %f (%f)\n', ticks/tic_end, tic_end);

%toc
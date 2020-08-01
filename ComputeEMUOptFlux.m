function [exitflag error flux cy_mt_ratio idv_opt] = ComputeEMUOptFlux(model, EMU, idv, EMU_met_known_mat, met_list_norm, WC_known_metabolites, initial_flux, initial_cy_mt_ratio)
global g_model g_EMU g_idv g_idv_known_arr g_idv_known_mat g_iteration_num g_met_list_norm g_WC_known_metabolites;
g_EMU = EMU;
g_idv = idv;
g_model = model;
g_met_list_norm = met_list_norm;
g_iteration_num = 0;
g_WC_known_metabolites = WC_known_metabolites;


g_idv_known_arr = EMU_met_known_mat(:);
g_idv_known_arr(isnan(g_idv_known_arr))=[];
g_idv_known_mat = EMU_met_known_mat;

[met_num, rxn_num] = size(model.S);

Aeq = model.S(model.met_extra == 0, :);
Aeq = [Aeq;model.equality_constraints];

Aeq = [Aeq zeros(size(Aeq,1),sum(sum(isnan(EMU_met_known_mat'))==0))];
Beq = zeros(size(Aeq,1),1);


A(5,8)=-1;A(5,9)=1;
b(5,1)=-0.001;

% A(5,8)=-1;A(5,9)=1;
% b(5,1)=-0.005;




% options = optimset('Algorithm','sqp', 'GradObj','on');
options = optimset('Algorithm','sqp', 'GradObj','off', 'MaxFunEvals', 5000);

%clc;
fprintf('Running fmincon...\n');
lb = [model.lb;zeros(size(initial_cy_mt_ratio))]; %flux lb and cy-mt ration of 0
ub = [model.ub;ones(size(initial_cy_mt_ratio))];  %flux ub and cy-mt ration of 1


[optimization_params, fval, exitflag] = fmincon(@opt_func, [initial_flux;initial_cy_mt_ratio],[],[], Aeq, Beq, lb, ub, [], options);
flux = optimization_params(1:g_model.rxn_num);
cy_mt_ratio = optimization_params(g_model.rxn_num+1:end);
fprintf('Finished fmincon\n');

fprintf('Iteration number = %d\n', g_iteration_num);


fprintf('Fluxes:\n');
for i=1:rxn_num
    fprintf('\tReaction %s (%d) = %.7f\n', model.rxns{i}, i, flux(i));
end
fprintf('CY_MT_Ratios:\n');
for i=1:length(WC_known_metabolites)
    x_cy = g_idv_known_mat(i,1);
    x_mt = g_idv_known_mat(i,2);
    if(((~isnan(x_cy)) && (~isnan(x_mt))))
       fprintf('\tMetabolite %s = %.7f\n', WC_known_metabolites{i}.met_name(1:end-2), cy_mt_ratio(i));
    end        
end
fprintf('\n\tRow constraints error = %f\n', sum(abs(Aeq*optimization_params)));

[idv_opt idv_d cycle_error] = ComputeEmuIDV(EMU, idv, g_idv_known_arr, flux);

fprintf('\nIDVs:\n');
for i=1:length(EMU.list)
    fprintf('\tEMU %s (%d) = %s\n', EMU.name{i}, i, DispIDV(idv_opt{i}));
end
e = ComputeError(g_WC_known_metabolites, idv_opt, g_idv_known_arr, g_idv_known_mat, g_EMU, g_model, cy_mt_ratio, flux);
error = e;
fprintf('\n\tIDV Error = %f\t cycle error %f\n', e, cycle_error+0);
fprintf('\n\t**************************  IDV Error = %f\t cycle error %f\n', e, cycle_error+0);



function [e] = opt_func(optimization_params)
global g_model g_EMU g_idv g_WC_known_metabolites g_idv_known_arr g_idv_known_mat g_iteration_num;
g_iteration_num = g_iteration_num+1;
flux = optimization_params(1:g_model.rxn_num);
cy_mt_ratio = optimization_params(g_model.rxn_num+1:end);
[idv idv_d cycle_error] = ComputeEmuIDV(g_EMU, g_idv, g_idv_known_arr, flux);
e = ComputeError(g_WC_known_metabolites, idv, g_idv_known_arr, g_idv_known_mat, g_EMU, g_model, cy_mt_ratio, flux);
% d = ComputeErrorDeriv(g_WC_known_metabolites, idv, idv_d, g_idv_known_arr, g_idv_known_mat, g_EMU, cy_mt_ratio, flux);

%  e = e+cycle_error
 
% fprintf('%f\n', e);



function e = ComputeError(g_WC_known_metabolites, idv, g_idv_known_arr, g_idv_known_mat, g_EMU, model, cy_mt_ratio, flux)
e = 0;
for i=1:size(g_idv_known_mat,1)    
    x_cy = g_idv_known_mat(i,1);
    x_mt = g_idv_known_mat(i,2);

    EMU_indices_cy = find(g_EMU.list(:,1)==x_cy);
    EMU_indices_mt = find(g_EMU.list(:,1)==x_mt);
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
       iterator_cy_mt_ratio = cy_mt_ratio(i);
    end    
    
    % compare known idv to the first EMU idv of this metabolite
    % it must be the first one from all EMU of the same metabolite, as the
    % first one contains all the carbons
    g_WC_known_metabolites{i}.idv_variance(g_WC_known_metabolites{i}.idv_variance<0.0001)=0.0001;
    try
    e = e + sum(((g_WC_known_metabolites{i}.idv(2:end)-((iterator_cy_mt_ratio*idv_cy(2:end))+((1-iterator_cy_mt_ratio)*idv_mt(2:end))) ).^2)./g_WC_known_metabolites{i}.idv_variance(2:end));
    catch
        sprintf('error on ComputeError');
    end
end
% e1=e;
% HF - updated from Won
% e = e + (((flux(1)-flux(2))-7.74)/0.016)^2;
% e = e + (((flux(16)-flux(17))-0.08)/0.02)^2;
% e = e + ((flux(27)-1.24)/0.03)^2;

% % LF - updated from Won
% e = e + (((flux(1)-flux(2))-1.68)/0.03)^2;
% e = e + (((flux(16)-flux(17))+0.84)/0.065)^2;
% e = e + ((flux(27)-0.18)/0.005)^2;

% add the distance from the measured net fluxes to the objective function
for(i=1:size(model.measured_net_fluxes_matrix,1))
    % for forward backward fluxes
    if(model.measured_net_fluxes_matrix(i,2)~=0)
        e = e + (((flux(model.measured_net_fluxes_matrix(i,1))-flux(model.measured_net_fluxes_matrix(i,2)))-model.measured_net_fluxes_matrix(i,3))/model.measured_net_fluxes_matrix(i,4))^2;
    else
        e = e + (((flux(model.measured_net_fluxes_matrix(i,1)))-model.measured_net_fluxes_matrix(i,3))/model.measured_net_fluxes_matrix(i,4))^2;
    end
end
    

function [c1 c2] = con_func(x)
    if (abs(x(1)) < 1e-2)
        c1 = 1;
    else
        c1 = 0;
    end
    
    c2 = [];
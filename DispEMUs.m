function DispEMUs(model, met_isotopomers, met_isotopomers_arr, v)

for i=v
    met = met_isotopomers_arr(i,1);
    met_isotope = met_isotopomers_arr(i,2);
   
    %fprintf('%d:\t', i);
    s = DispEMU(model.mets{met}, met_isotopomers{met}(:, met_isotope));    
    fprintf('%s,\t',s);
end
fprintf('\n');

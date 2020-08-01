function idv = CreateIDV(EMU, extra_met_isotopomers)

n = length(EMU.list);
idv = cell(n,1);

for i=1:length(EMU.list)
    met = EMU.list(i,1);
    met_emu_ind = EMU.list(i,2);
    if ~isempty(extra_met_isotopomers{met})
        t = EMU.met{met}(:, met_emu_ind);
        v = extra_met_isotopomers{met}.atom_idv(t==1);
        enrichment = extra_met_isotopomers{met}.enrichment;
        a = sum(v);
        
        idv{i} = zeros(1, EMU.size(i)+1);
        
        if (a>0)        
            idv{i}(a+1) = 1*enrichment;
            idv{i}(1) = 1*(1-enrichment);
        else
            idv{i}(1) = 1;            
        end              
    end
end

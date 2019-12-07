function s = DispEMU(met_name, emu)
v = find(emu);
s = '';
for i=v'
    s = sprintf('%s%d', s, i);
end
s = sprintf('%s_%s', met_name, s);
    
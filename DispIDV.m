function s = DispIDV(v)
s = '[';
for i=1:length(v)
    if (i>1)
        s = sprintf('%s, ', s);
    end
    s = sprintf('%s%.3f', s, v(i)+0);
end
s = sprintf('%s]',s);
    
function [ output_args ] = bar_label( D, L )

    hBar = bar(D, 'stacked');
    
    yd = get(hBar, 'YData');
    
    %set(gca, 'XTick', [-1:1:size(D,1)]); %%?
    %xt = get(gca, 'XTick'); 
    xt = [1:size(D,1)];
    
    barbase = cumsum([zeros(size(D,1),1) D(:,1:end-1)],2);
    joblblpos = D/2 + barbase;
    
    for k1 = 1:size(D,1)
        v = L(k1,:);
        yjob = cell(0);
        for i=1:length(v)
            yjob{i} = round(v(i)*100)/100;
            if (v(i) < 0.05)
                yjob{i} = '';
            end
                
        end
    text(xt(k1)*ones(1,size(D,2)), joblblpos(k1,:), yjob, 'HorizontalAlignment','center', 'Color','w', 'FontSize', 16)
    end
end


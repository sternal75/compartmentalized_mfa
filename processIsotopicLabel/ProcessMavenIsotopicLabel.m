function [met_list_norm, met_name_arr, sample_list_short] = ProcessMavenIsotopicLabel(file_name, token_number)
[d,t] = xlsread(file_name);

met_id = d(:,2);
met_mass = d(:,4);
met_RT = d(:,5);
met_list = t(2:end,9);
met_note = t(2:end,8);
met_mat = d(:, 13:end);
sample_list = t(1,14:end);

% Remove () from metabolite name
for i=1:length(met_list)
    s = met_list{i};
    n = strfind(s, '(');
    if isempty(n)
        continue;
    end
    n = n(end);
    s = s(1:n-2);
    
    % remove some characters
    s = strrep(s, '?', '');
    met_list{i} = s;
end

% Extract sample ID 
sample_list2 = cell(0);
for i=1:length(sample_list)
    s = sample_list{i};
    a=strfind(s, '_');
    
    sample_list2{i} = s;
    if length(a)>=2
        sample_list2{i} = s((a(2)+1):end);
    end
end

% Remove blans from notes 
for i=1:length(met_note)
    if length(met_note{i}) < 1
        met_note{i} = 'C12 PARENT';
    end
end


% Remove blanks
v = ones(length(sample_list2),1);
n = strmatch('Blank', sample_list2);
v(n) = 0;
sample_list2 = sample_list2(v==1);
met_mat = met_mat(:, v==1);

% Find sample replicates
[a,b,c]=unique(sample_list2);
% a
if (1==1)
    curr = c(1);
    mat = [];
    ind = 1;
    sample_mat = [];
    sample_mat(ind, 1) = 1;
    sample_list_short{1} = a{curr};
    for i=2:length(c)
        if (curr ~= c(i))
            ind = ind + 1;
            sample_list_short{ind} = a{c(i)};
        end
        curr = c(i);
        sample_mat(ind, i) = 1;
    end
else
    n = length(a);
    v = [];
    for i=1:n
        t = find(c==i);
        sample_mat(i, t) = 1;
        v(i) = str2num(a{i});
    end
    [a,b] = sort(v);
    sample_mat = sample_mat(b, :);
    
    for i=1:length(a)
        sample_list_short{i} = num2str(a(i));
    end
    
end

% Iterate over metabolites and normalize data
v = strmatch('C12', met_note);
met_list_norm = cell(0);
met_name_arr = cell(0);
for i=1:length(v)
    b = v(i);
    if i<length(v)
        e = v(i+1)-1;
    else
        e = size(met_mat,1);
    end
        
    m = met_mat(b:e,:);
    t = sum(m,1);
    m2 = t;
    met_list_norm{i}.median_intensity = median(t);
    m = m ./ repmat(t, size(m,1),1);
%   Alon - added the below line as in somecases when normalizing 
%   the metabolites, the values were devided by 0 and become NAN. In that
%   case it should be zero
    m(isnan(m)) = 0 ;
    
    met_list_norm{i}.met_name = met_list{b};
    if (strcmp(met_list_norm{i}.met_name, 'N-Acetylputrescine')==1)
        a=a;
    end
    met_name_arr{i} = met_list{b};
    met_list_norm{i}.data = zeros(size(m,1), size(sample_mat,1));
    met_list_norm{i}.data_total = zeros(1, size(sample_mat,1));
    met_list_norm{i}.data_total_m = cell(0);
    for j=1:size(sample_mat,1)
        u = find(sample_mat(j,:));
        
        met_list_norm{i}.data(:, j) = mean(m(:, u)')';
        met_list_norm{i}.var(:, j)  = var(m(:, u)')';        
        met_list_norm{i}.data_total(j) = mean(m2(u)')';
        met_list_norm{i}.data_total_m{j} = m2(u);
    end
    met_list_norm{i}.data_total = met_list_norm{i}.data_total ./ max(met_list_norm{i}.data_total);
end

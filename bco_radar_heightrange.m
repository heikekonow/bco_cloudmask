clear
path = '/pool/OBS/BARBADOS_CLOUD_OBSERVATORY/Level_1/B_Reflectivity/Version_2/';
radarname = {'KATRIN', 'MBR'};

for i=1:length(radarname)
    files = listFiles([path '*' radarname{i} '*.nc']);
    a{i} = cellfun(@(x) regexp(x, '__'), files, 'uni', false);
    heightrange{i} = cellfun(@(x,y) x(y(4)+2:y(5)-1), files, a{i}, 'uni', false);
    unique_height{i} = unique(heightrange{i});

    for j=1:length(unique_height{i})
        radarfiles{i,j} = listFiles([path '*' radarname{i} '*' unique_height{i}{j} '*.nc']);
        dates{i,j} = cell2mat(cellfun(@(x) x(a{i}{j}(5)+2:end-3), radarfiles{i,j}, 'uni', false));
    end
end

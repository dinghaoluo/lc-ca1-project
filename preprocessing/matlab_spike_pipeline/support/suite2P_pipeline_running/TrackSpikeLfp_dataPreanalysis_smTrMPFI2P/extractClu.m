function Clu = extractClu(data_2p)
% extract the cells that are manually curated

    % indice of manually selected cells
    indSel = find(data_2p.iscell(:,1) == 1);
    
    % total number of cells
    totNoCells = size(data_2p.iscell,1);
    
    % label all the cells that are merged into one cell
    merged = zeros(1,totNoCells);
    for i = 1:length(indSel)
        if(isfield(data_2p.stat{indSel(i)},'imerge') && ~isempty(data_2p.stat{indSel(i)}.imerge))
            merged(data_2p.stat{indSel(i)}.imerge+1) = 1;
        end
    end
    indMerged = find(merged == 1);
    
    % final cell list
    indSel = setdiff(indSel, indMerged);
    totNoCells = length(indSel);
    
    % getting the statistics of the cells
    Clu = struct(...
                'localClu', zeros(1,totNoCells),...
                'ypix', {cell(1,totNoCells)},...
                'xpix', {cell(1,totNoCells)},...
                'lam', {cell(1,totNoCells)},...
                'med', zeros(2,totNoCells),...
                'footprint', zeros(1,totNoCells),...
                'mrs', zeros(1,totNoCells),...
                'compact', zeros(1,totNoCells),...
                'solidity', zeros(1,totNoCells),...
                'npix', zeros(1,totNoCells),...
                'radius', zeros(1,totNoCells),...
                'aspect_ratio', zeros(1,totNoCells),...
                'npix_norm', zeros(1,totNoCells),...
                'skew', zeros(1,totNoCells),...
                'std', zeros(1,totNoCells),...
                'chan2_prob', zeros(1,totNoCells));
            
    for i = 1:totNoCells
        Clu.localClu(i) = indSel(i);
        Clu.ypix{i} = data_2p.stat{indSel(i)}.ypix;
        Clu.xpix{i} = data_2p.stat{indSel(i)}.xpix;
        Clu.lam{i} = data_2p.stat{indSel(i)}.lam;
        Clu.med(:,i) = data_2p.stat{indSel(i)}.med';
        Clu.footprint(i) = data_2p.stat{indSel(i)}.footprint;
        Clu.mrs(i) = data_2p.stat{i}.mrs;
        Clu.compact(i) = data_2p.stat{indSel(i)}.compact;
        Clu.solidity(i) = data_2p.stat{indSel(i)}.solidity;
        Clu.npix(i) = data_2p.stat{indSel(i)}.npix;
        Clu.radius(i) = data_2p.stat{indSel(i)}.radius;
        Clu.aspect_ratio(i) = data_2p.stat{indSel(i)}.aspect_ratio;
        Clu.npix_norm(i) = data_2p.stat{indSel(i)}.npix_norm;
        Clu.skew(i) = data_2p.stat{indSel(i)}.skew;
        Clu.std(i) = data_2p.stat{indSel(i)}.std;
        Clu.chan2_prob(i) = data_2p.stat{indSel(i)}.chan2_prob;
    end
end

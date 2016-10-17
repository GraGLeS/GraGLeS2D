function [ebsd_cor,grains] =  correctEBSD(ebsd,opts)


%% detect grains based on segangle
grains = calcGrains(ebsd,'angle',opts.segAngle*degree,'augmentation','convexhull','silent');    
   
%sort out small grains = at least nFilter measurements per grain
if opts.nFilter > 0
    
    for iter = 1:opts.iter
        % discriminate on basis of nFilter
        fprintf('\nFiltering with at least %d points per grain\n',opts.nFilter);
        grains = grains(grainSize(grains) >= opts.nFilter);
        ebsd_cor = get(grains,'EBSD');

        % exctrapolate nearest neighbour
        ebsd_cor = fill(ebsd_cor,extend(ebsd),opts.stepsize);

        % calculate grains from filtered ebsd
        fprintf('Recalculating grains...');
        grains = calcGrains(ebsd_cor,'angle',opts.segAngle*degree); %'convexhull','augmentation','silent'
    end
else
    % do not apply filtering
    ebsd_cor = ebsd;
end

%% smooth grains
    
if opts.sm > 0
    grains = smooth(grains,opts.sm);
end

 fprintf('done (%d grains)\n',numel(grains));

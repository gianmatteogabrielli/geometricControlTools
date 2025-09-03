function jumpTimes = generateClusteredJumps(t_end, mainJumps, jumpsPerCluster, clusterSpread)
    rng(1)
% generateClusteredJumps genera istanti di salto raggruppati attorno ai mainJumps
%
% INPUTS:
%   t_end           - tempo finale
%   mainJumps       - array con istanti principali di salto
%   jumpsPerCluster - numero di salti da generare attorno a ciascun mainJump (default: 5)
%   clusterSpread   - deviazione standard della distribuzione attorno a ciascun mainJump (default: 0.5)
%
% OUTPUT:
%   jumpTimes       - array contenente tutti gli istanti di salto generati

    if nargin < 3
        jumpsPerCluster = 5;
    end
    if nargin < 4
        clusterSpread = 0.5; % secondi
    end

    jumpTimes = [];

    for i = 1:length(mainJumps)
        mu = mainJumps(i);
        % Genera tempi gaussiani attorno a mu, con spread
        cluster = mu + clusterSpread * randn(1, jumpsPerCluster);
        % Mantieni solo quelli nel range [0, t_end]
        cluster = cluster(cluster >= 0 & cluster <= t_end);
        jumpTimes = [jumpTimes, cluster];
    end

    % Ordina i tempi di salto
    jumpTimes = sort(jumpTimes);
end

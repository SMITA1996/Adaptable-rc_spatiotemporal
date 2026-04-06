%% ==================== PQ MODEL: STORE SNAPSHOTS IN WORKSPACE ====================
clear; clc;

% ---------------- Parameters ----------------
N = 40;                     % lattice size
q = 0.92;                    % positive feedback strength
p_values = 0.01:0.001:1;%bifurcation
T = 2e6;                     % number of update attempts per p

Results = struct();         

figure('Color','w'); colormap(gray);

%% ---------------- Simulation Loop ----------------
for k = 1:length(p_values)
    p = p_values(k);
    
    % Initialize lattice with initial density 0.5
    A = rand(N,N) < 0.5;

    % --- Monte Carlo Dynamics ---
    for step = 1:T
        i = randi(N);
        j = randi(N);

        if A(i,j) == 1
            [in,jn] = random_neighbor(i,j,N);

            if A(in,jn) == 0
                if rand < p
                    A(in,jn) = 1;
                else
                    A(i,j) = 0;
                end
            else
                if rand < q
                    idx = neighbor_of_pair(i,j,in,jn,N);
                    A(idx) = 1;
                elseif rand < (1-p)
                    A(i,j) = 0;
                end
            end
        end
    end

    % Store result in workspace
    Results(k).p = p;
    Results(k).A = A;
    Results(k).mean_density = mean(A(:));

    % Display snapshot 
    subplot(1,length(p_values),k);
    imagesc(A); axis square off;
    
end
function [in,jn] = random_neighbor(i,j,N)
dirs = [-1 0;1 0;0 -1;0 1];
d = dirs(randi(4),:);
in = mod(i-1+d(1),N)+1;
jn = mod(j-1+d(2),N)+1;
end

function idx = neighbor_of_pair(i,j,in,jn,N)
neigh = [
    i-1 j; i+1 j; i j-1; i j+1; 
    in-1 jn; in+1 jn; in jn-1; in jn+1
];
neigh = mod(neigh-1,N)+1;        % periodic BC
neigh = unique(neigh,'rows');    % remove duplicates
neigh = neigh(~(neigh(:,1)==i & neigh(:,2)==j),:);
neigh = neigh(~(neigh(:,1)==in & neigh(:,2)==jn),:);

r = randi(size(neigh,1));
idx = sub2ind([N N],neigh(r,1),neigh(r,2));
end


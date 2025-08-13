function h = daboxplot(Y,varargin)
% daboxplot - version modifiée pour tracer des lignes entre 3 conditions dans deux groupes distincts.
% ...

h = struct;
p = inputParser;

% Spécification des options par défaut
addOptional(p, 'groups', []);
addOptional(p, 'fill', 1);
addOptional(p, 'colors', get(gca,'colororder'));
addOptional(p, 'whiskers', 1);
addOptional(p, 'scatter', 0);
addOptional(p, 'scattersize', 15);
addOptional(p, 'scattercolors', {'k','w'});
addOptional(p, 'flipcolors', 0);
addOptional(p, 'scatteralpha', 1);
addOptional(p, 'jitter', 1);
addOptional(p, 'mean', 0);
addOptional(p, 'outliers', 1);
addOptional(p, 'outfactor', 1.5);
addOptional(p, 'outsymbol', 'rx');
addOptional(p, 'boxalpha', 1);
addOptional(p, 'boxspacing', 1);
addOptional(p, 'boxwidth', 1);
addOptional(p, 'linkline', 0);
addOptional(p, 'withinlines', 0);
addOptional(p, 'linklinesGroups', 0); % Nouveau paramètre
addOptional(p, 'xtlabels', []);
addOptional(p, 'legend', []);

% Analyse des paramètres
parse(p, varargin{:});
confs = p.Results;

% Détection du nombre de groupes et conditions
if iscell(Y)
    num_groups = numel(Y);
    Y = cell2mat(cellfun(@(x) reshape(x, [], size(x, 2)), Y, 'UniformOutput', false));
    Gi = repelem(1:num_groups, size(Y, 1) / num_groups)';
else
    if isempty(confs.groups)
        Gi = ones(size(Y, 1), 1);
        num_groups = 1;
    else
        [Gi, ~] = grp2idx(confs.groups);
        num_groups = numel(unique(Gi));
    end
end

num_conditions = size(Y, 2);
cpos = 1:num_conditions;

% Calcul des positions des boîtes
if num_groups == 1
    gpos = cpos;
    box_width = 1 / 3 * confs.boxwidth;
else
    box_width = (2 / 3) / (num_groups + 1) * confs.boxwidth;
    loc_sp = (box_width / 3) * confs.boxspacing;
    gpos = arrayfun(@(x) cpos + (x - (num_groups + 1) / 2) * (box_width + loc_sp), 1:num_groups, 'UniformOutput', false);
    gpos = cell2mat(gpos(:));
end

h.gpos = gpos;
h.cpos = cpos;

% Création des boîtes et des lignes par groupes
for g = 1:num_groups
    group_data = Y(Gi == g, :);
    pt = prctile(group_data, [25 50 75]);
    means = mean(group_data, 1);

    % Positions des boîtes pour chaque condition
    g_indices = (1:num_conditions) + (g - 1) * num_conditions;

    % Tracer les boîtes
    for k = 1:num_conditions
        box_pos = gpos(g_indices(k));
        box_x = [box_pos - box_width / 2, box_pos + box_width / 2];
        fill(box_x([1 2 2 1]), pt([1 1 3 3], k), confs.colors(g, :), 'FaceAlpha', confs.boxalpha, 'EdgeColor', 'none');
        hold on;
        plot(box_x, pt([2 2], k), 'k-', 'LineWidth', 1.5);
    end

    % Tracer les lignes reliant les points des conditions pour chaque participant
    if confs.linklinesGroups
        for i = 1:size(group_data, 1)
            plot(gpos(g_indices), group_data(i, :), '-', 'Color', [0.7 0.7 0.7]);
        end
    end
end

% Gestion des étiquettes et des limites
if ~isempty(confs.xtlabels)
    set(gca, 'XTick', cpos, 'XTickLabels', confs.xtlabels);
else
    set(gca, 'XTick', cpos);
end
xlim([gpos(1) - box_width, gpos(end) + box_width]);

% Légende
if ~isempty(confs.legend)
    legend(confs.legend);
end

end

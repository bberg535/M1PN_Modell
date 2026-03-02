function registry = mm_model_registry(cfg)
%MM_MODEL_REGISTRY Enumerate models and default orders from config.

registry = struct('name', {}, 'order', {}, 'nMom', {}, 'label', {});
idx = 0;
for i = 1:numel(cfg.models.families)
    fam = cfg.models.families{i};
    if ~isfield(cfg.models.orders, fam)
        continue;
    end
    orders = cfg.models.orders.(fam);
    for j = 1:numel(orders)
        ord = orders(j);
        m = mm_build_model(fam, ord, cfg.models);
        idx = idx + 1;
        registry(idx).name = fam;
        registry(idx).order = ord;
        registry(idx).nMom = m.nMom;
        registry(idx).label = sprintf('%s-%d', fam, ord);
    end
end

end

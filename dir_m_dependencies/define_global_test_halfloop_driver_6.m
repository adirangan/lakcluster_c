global_prefix_method_ = { ...
 'tsne00_isosplit5' ...
,'tsne50_isosplit5' ...
,'umap00_default' ...
,'umap00_isosplit5' ...
,'umap00_hdbscan' ...
,'louvain00_default' ...
,'hnbtZRgumb' ...
,'spectral_isosplit5' ...
,'tsne00pr_isosplit5' ...
,'tsne50pr_isosplit5' ...
,'umap00pr_default' ...
,'umap00pr_isosplit5' ...
,'umap00pr_hdbscan' ...
,'louvain00pr_default' ...
,'hnbrtZRgumb' ...
,'hnbr0tZRgumb' ...
};
global_n_method = numel(global_prefix_method_);
global_index_method_proj_off_ = ...
[ ...
 efind(cellfun(@(x)strcmp(x,'tsne00_isosplit5'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'tsne50_isosplit5'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'umap00_default'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'umap00_isosplit5'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'umap00_hdbscan'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'louvain00_default'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'hnbtZRgumb'),global_prefix_method_)) ...
];
global_index_method_proj_0on_ = ...
[ ...
 efind(cellfun(@(x)strcmp(x,'spectral_isosplit5'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'tsne00pr_isosplit5'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'tsne50pr_isosplit5'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'umap00pr_default'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'umap00pr_isosplit5'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'umap00pr_hdbscan'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'louvain00pr_default'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'hnbrtZRgumb'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'hnbr0tZRgumb'),global_prefix_method_)) ...
];
global_index_method_sub_ = ...
[ ...
 efind(cellfun(@(x)strcmp(x,'umap00_default'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'louvain00_default'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'hnbtZRgumb'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'umap00pr_default'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'louvain00pr_default'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'hnbrtZRgumb'),global_prefix_method_)) ...
,efind(cellfun(@(x)strcmp(x,'hnbr0tZRgumb'),global_prefix_method_)) ...
];
%%%%%%%%;
global_legend_method_ = { ...
 'tsne00_isosplit5' ...
,'tsne50_isosplit5' ...
,'umap00_default' ...
,'umap00_isosplit5' ...
,'umap00_hdbscan' ...
,'louvain00_default' ...
,'halflooprc_c' ...
,'spectral_isosplit5' ...
,'tsne00pr_isosplit5' ...
,'tsne50pr_isosplit5' ...
,'umap00pr_default' ...
,'umap00pr_isosplit5' ...
,'umap00pr_hdbscan' ...
,'louvain00pr_default' ...
,'halflooprcpr_c' ...
,'halfloopr0pr_c' ...
};
for nl=0:numel(global_legend_method_)-1; global_legend_method_{1+nl}(strfind(global_legend_method_{1+nl},'_')) = ' '; end;
global_nrank_method_ = cell(global_n_method,1);
for nmethod=0:global_n_method-1;
global_nrank_method_{1+nmethod} = [];
if ~isempty(intersect(nmethod,global_index_method_proj_0on_));
global_nrank_method_{1+nmethod} = 1;
end;%if ~isempty(intersect(nmethod,global_index_method_proj_0on_));
end;%for nmethod=0:global_n_method-1;
global_symbol_method_ = { ...
 'co-' ...
,'bo-' ...
,'kx-' ...
,'ko-' ...
,'k^-' ...
,'rs-' ...
,'mh-' ...
,'go-' ...
,'co-' ...
,'bo-' ...
,'kx-' ...
,'ko-' ...
,'k^-' ...
,'rs-' ...
,'mh-' ...
,'mp-' ...
};
%%%%%%%%;
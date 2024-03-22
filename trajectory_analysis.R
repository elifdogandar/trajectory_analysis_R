
# Trajectory analysis with Monocle3 with a Seurat object

library(Seurat)
library(monocle3)
library(tidyverse)
library("biomaRt")
library(SeuratDisk)
library(dplyr)
library(patchwork)

seurat_obj <-readRDS(file = "seurat_obj.rds")

DimPlot(seurat_obj, reduction = "umap",
        label = TRUE, group.by = "seurat_clusters",
        label.size = 6, pt.size =2.5)

# Subsetting to some clusters only

subsetted <- subset(x = seurat_obj,
                    subset = seurat_clusters %in% c('1','7','0','10','11') )

DimPlot(subsetted, reduction = "umap",
        label = TRUE, group.by = "seurat_clusters",
        label.size = 6, pt.size =2.5)+
  labs(title = "Subsetted cells")

FeaturePlot(subsetted, features = "X_score")

# Transforming from Seurat to cell data set
cds <- SeuratWrappers::as.cell_data_set(subsetted)
cds <- estimate_size_factors(cds)

# For plot_cells to work with gene names
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

cds <- cluster_cells(cds,resolution = 1/200) ## Step 4: Cluster the cells

# Show the clusters calculated by monocle3
plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)

# Show partitions calculated by monocle3
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)

cds <- learn_graph(cds) ## Step 5: Learn a graph
cds <- order_cells(cds) ## Step 6: Order cells 

# Plotting the trajectory graph white circles are leaves dark circles are forks
plot_cells(cds,
           show_trajectory_graph = T,
           color_cells_by = "seurat_clusters",
           label_cell_groups = FALSE,
           group_label_size = 6,
           label_groups_by_cluster = F,
           cell_size = 0.8,
           graph_label_size = 3)

# Plotting calculated pseudotime
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 0.7)+
  labs(title = "resolution = 1/200")

# Saving pseudotime as metadata
cds$monocle3_pseudotime <- pseudotime(cds)

# Boxplots of pseudotime for each cluster
data.pseudo1 <- as.data.frame(colData(cds))

ggplot(data.pseudo1, aes(monocle3_pseudotime, seurat_clusters,
                         fill = seurat_clusters)) + geom_boxplot()

ggplot(data.pseudo1, aes(monocle3_pseudotime,
                         reorder(seurat_clusters, monocle3_pseudotime),
                         fill = seurat_clusters)) + geom_boxplot() 
+labs(y='resolution = 1/200')

ggsave("pseudotime_by_cluster.png", dpi = 600)




# Regression analysis
gene_fits <- fit_models(cds, model_formula_str = "~Group")
fit_coefs <- coefficient_table(gene_fits)
Group_terms <- fit_coefs %>% filter(term == "Group_WT")
Group_terms %>% filter (q_value < 0.05) %>%
  dplyr::select(gene_id, term, q_value, estimate)


genes <- c("Rpl39",
           "Tmsb4x",
           "Dnaja1")

cds_subset <- cds[rowData(cds)$gene_short_name %in% genes,]
plot_genes_violin(cds_subset, group_cells_by="seurat_clusters", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# Should we add batch effect to the model too?
group_batch_models <- fit_models(cds_subset,
                                model_formula_str = "~Group + batch",
                                expression_family="negbinomial")
group_models <- fit_models(cds_subset,
                          model_formula_str = "~Group",
                          expression_family="negbinomial")
compare_models(group_batch_models, group_models) %>% select(gene_short_name, q_value)


# Differential gene expression throughout the trajectory

cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=8)

cds_pr_test_res_ok <- arrange(cds_pr_test_res, q_value)%>% filter(status == "OK")
cds_pr_test_res_ok%>% head()
cds_pr_test_res2 <- arrange(cds_pr_test_res, q_value)%>% filter(status == "OK")%>% top_n(-12)


top12 <- row.names( cds_pr_test_res2)[1:12]

# Plotting top 12 trajectory dependent genes
plot_cells(cds, genes=top12,
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           cell_size = 0.4)
RidgePlot(subsetted, features = top12[1:4],
          sort = T,ncol = 2)

# Finding modules of co-regulated genes
pr_deg_ids <- row.names(subset(cds_pr_test_res_ok, q_value < 0.005))

cds <- preprocess_cds(cds) # necessary for gene modules
gene_modules <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))
table(gene_modules$module)


cell_groups <- data.frame(cell = row.names(colData(cds)),
                          cell_group = colData(cds)$seurat_clusters)
agg_mat <- aggregate_gene_expression(cds,
                                     gene_group_df = gene_modules,
                                     cell_group_df = cell_groups)
dim(agg_mat)

row.names(agg_mat) <- paste0("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column",
                   treeheight_row = 0,
                   treeheight_col = 0,
                   clustering_method="ward.D2")

gm <- gene_modules[which(gene_modules$module %in% c(1, 1)),]
plot_cells(cds,
           genes=gm$id[1:15],
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           trajectory_graph_color = "grey60")

cds_subset <-  cds[rowData(cds)$gene_short_name %in% top12,]
cds_subset <- order_cells(cds_subset)

plot_genes_in_pseudotime(cds_subset,
                         color_cells_by="seurat_clusters",
                         min_expr=0.5)

plot_genes_violin(cds_subset, group_cells_by="seurat_clusters", ncol=4)

plot_cells(cds, genes='Hspa1b',
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           cell_size = 0.4)


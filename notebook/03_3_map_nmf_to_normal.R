# run harmony if needed
#' R 4.3.3

library(Seurat)
library(tidyverse)
library(here)
library(schard)
options(future.globals.maxSize = 100 * 1024^3)

tumor <- h5ad2seurat(here("output/9_annotated_subcluster_embedding_tumor_nmf_added.h5ad"))
tumor$NMF <- apply(tumor@meta.data[, paste0("NMF_Module_", 1:6)], 1, function(x) {
    paste0("NMF_Module", which.max(x))
})


tumor <- tumor |>
    NormalizeData(verbose = FALSE) |>
    FindVariableFeatures(verbose = FALSE) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(npcs = 20, verbose = FALSE)
# DimPlot(tumor, group.by = c("cell_type", "sample_id"))


# # fetal -------------------------#
# adata_f <- h5ad2seurat(here("public/allen_fetal_nmf_filtered.h5ad"))
# adata_f <- adata_f |>
#     NormalizeData(verbose = FALSE) |>
#     FindVariableFeatures(verbose = FALSE) |>
#     ScaleData(verbose = FALSE) |>
#     RunPCA(npcs = 20, verbose = FALSE) |>
#     RunUMAP(dims = 1:20, reduction = "pca", return.model = TRUE, verbose = FALSE)
# DimPlot(adata_f, reduction = "umap", group.by = "cell_type")

# ## query
# anchors <- FindTransferAnchors(reference = adata_f, query = tumor, dims = 1:20, reference.reduction = "pca", k.filter = NA)

# predictions <- TransferData(anchorset = anchors, refdata = adata_f$cell_type, dims = 1:20)

# tumor <- AddMetaData(tumor, metadata = predictions)
# tumor <- MapQuery(
#     anchorset = anchors, reference = adata_f, query = tumor,
#     refdata = list(celltype = "cell_type"), reference.reduction = "pca", reduction.model = "umap"
# )

# ## plot
# p1 <- DimPlot(adata_f,
#     reduction = "umap", group.by = "cell_type",
#     label = TRUE, label.size = 3, repel = TRUE, alpha = 0.1,
#     cols = rep("grey", length(unique(adata_f$cell_type)))
# ) +
#     NoLegend()
# p2 <- DimPlot(tumor,
#     reduction = "ref.umap", group.by = "predicted.id",
#     label = TRUE, label.size = 3, repel = TRUE
# ) +
#     NoLegend()
# combined_plot <- p1 +
#     geom_point(
#         data = p2$data,
#         aes(x = refUMAP_1, y = refUMAP_2, color = predicted.id),
#         alpha = 0.8, size = 2
#     ) +
#     scale_color_discrete()
# ggsave(here("output/fetal_proj.png"), plot = combined_plot, width = 10, height = 5)





# # adult -------------------------#
# adata_a <- h5ad2seurat(here("public/allen_adult_nmf_filtered.h5ad"))
# adata_a <- adata_a |>
#     NormalizeData(verbose = FALSE) |>
#     FindVariableFeatures(verbose = FALSE) |>
#     ScaleData(verbose = FALSE) |>
#     RunPCA(npcs = 20, verbose = FALSE) |>
#     RunUMAP(dims = 1:20, reduction = "pca", return.model = TRUE, verbose = FALSE)
# DimPlot(adata_a, reduction = "umap", group.by = "cell_type")


# anchors <- FindTransferAnchors(reference = adata_a, query = tumor, dims = 1:20, reference.reduction = "pca", k.filter = NA)

# predictions <- TransferData(anchorset = anchors, refdata = adata_a$cell_type, dims = 1:20)

# tumor <- AddMetaData(tumor, metadata = predictions)
# tumor <- MapQuery(
#     anchorset = anchors, reference = adata_a, query = tumor,
#     refdata = list(celltype = "cell_type"), reference.reduction = "pca", reduction.model = "umap"
# )

# ## plot
# p1 <- DimPlot(adata_a,
#     reduction = "umap", group.by = "cell_type",
#     label = TRUE, label.size = 3, repel = TRUE, alpha = 0.1,
#     cols = rep("grey", length(unique(adata_a$cell_type)))
# ) +
#     NoLegend()
# p2 <- DimPlot(tumor,
#     reduction = "ref.umap", group.by = "predicted.id",
#     label = TRUE, label.size = 3, repel = TRUE
# ) +
#     NoLegend()
# combined_plot <- p1 +
#     geom_point(
#         data = p2$data,
#         aes(x = refUMAP_1, y = refUMAP_2, color = predicted.id),
#         alpha = 0.8, size = 2
#     ) +
#     labs(title = "Normal adult mouse projection") +
#     scale_color_discrete()
# ggsave(here("output/adult_proj.png"), plot = combined_plot, width = 10, height = 5)

# ## plot NMF
# p1 <- DimPlot(adata_a,
#     reduction = "umap", group.by = "cell_type",
#     label = TRUE, label.size = 3, repel = TRUE, alpha = 0.1,
#     cols = rep("grey", length(unique(adata_a$cell_type)))
# ) +
#     NoLegend()
# p2 <- DimPlot(tumor,
#     reduction = "ref.umap", group.by = "NMF",
#     label = TRUE, label.size = 3, repel = TRUE
# ) +
#     NoLegend()
# combined_plot <- p1 +
#     geom_point(
#         data = p2$data,
#         aes(x = refUMAP_1, y = refUMAP_2, color = NMF),
#         alpha = 0.8, size = 2
#     ) +
#     scale_color_discrete() +
#     labs(title = "Normal adult mouse projection") +
#     theme(legend.position = "bottom")
# ggsave(here("output/adult_proj_nmf.png"), plot = combined_plot, width = 10, height = 5)

# plot_data <- as.data.frame.matrix(table(tumor$predicted.id, tumor$NMF)) |>
#     tibble::rownames_to_column("predicted_id") |>
#     tidyr::pivot_longer(-predicted_id,
#         names_to = "NMF",
#         values_to = "count"
#     )

# ggplot(plot_data, aes(x = NMF, y = predicted_id, fill = count)) +
#     geom_tile() +
#     scale_fill_viridis_c() +
#     theme_minimal() +
#     labs(x = "NMF Module", y = "Predicted ID")

# integrate ----------------------------------- #

adata_a <- h5ad2seurat(here("public/allen_adult_nmf_filtered.h5ad"))
adata_a$merged_cell_type <- paste0("A-", adata_a$cell_type)
adata_a <- adata_a |>
    NormalizeData(verbose = FALSE) |>
    FindVariableFeatures(verbose = FALSE) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(npcs = 20, verbose = FALSE) |>
    RunUMAP(dims = 1:20, reduction = "pca", return.model = TRUE, verbose = FALSE)

adata_f <- h5ad2seurat(here("public/allen_fetal_nmf_filtered.h5ad"))
adata_f$merged_cell_type <- paste0("F-", adata_f$cell_type)
adata_f <- adata_f |>
    NormalizeData(verbose = FALSE) |>
    FindVariableFeatures(verbose = FALSE) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(npcs = 20, verbose = FALSE) |>
    RunUMAP(dims = 1:20, reduction = "pca", return.model = TRUE, verbose = FALSE)



features <- SelectIntegrationFeatures(list(adata_a, adata_f))
anchors <- FindIntegrationAnchors(object.list = list(adata_a, adata_f), anchor.features = features)
adata <- IntegrateData(anchorset = anchors)
DefaultAssay(adata) <- "integrated"

adata <- adata |>
    FindVariableFeatures() |>
    ScaleData() |>
    RunPCA(npcs = 20) |>
    RunUMAP(reduction = "pca", dims = 1:20, return.model = TRUE)

anchors <- FindTransferAnchors(reference = adata, query = tumor, dims = 1:20, reference.reduction = "pca", k.filter = NA)

predictions <- TransferData(anchorset = anchors, refdata = adata$cell_type, dims = 1:20)

tumor <- AddMetaData(tumor, metadata = predictions)
tumor <- MapQuery(
    anchorset = anchors, reference = adata, query = tumor,
    refdata = list(celltype = "merged_cell_type"), reference.reduction = "pca", reduction.model = "umap"
)

Sys.getpid()
saveRDS(adata, here("public/integrated_ref.rds"))
saveRDS(tumor, here("public/tumor_mapped.rds"))
save.image(here("public/workspace.RData"))
getwd()

# cell_type or merged_cell_type
# Plot
p1 <- DimPlot(adata,
    reduction = "umap", group.by = "cell_type",
    label = TRUE, label.size = 3, repel = TRUE, alpha = 0.1,
    cols = rep("grey", length(unique(adata$cell_type)))
) +
    NoLegend()
p2 <- DimPlot(tumor,
    reduction = "ref.umap", label.size = 3, repel = TRUE
) +
    NoLegend()

combined_plot <- p1 +
    geom_point(
        data = p2$data,
        aes(x = refUMAP_1, y = refUMAP_2),
        alpha = 0.4, size = 1, , color = "red"
    ) +
    labs(title = "Normal adult&fetal mouse projection") +
    theme(legend.position = "bottom")
ggsave(here("output/integ_proj.png"), plot = combined_plot, width = 10, height = 5)


## plot NMF
tumor[[]] |> View()
p1 <- DimPlot(adata,
    reduction = "umap", group.by = "cell_type",
    label = TRUE, label.size = 3, repel = TRUE, alpha = 0.1,
    cols = rep("grey", length(unique(adata$cell_type)))
) +
    NoLegend()
p2 <- DimPlot(tumor,
    reduction = "ref.umap", group.by = "NMF",
    label = TRUE, label.size = 3, repel = TRUE
) +
    NoLegend()
combined_plot <- p1 +
    geom_point(
        data = p2$data,
        aes(x = refUMAP_1, y = refUMAP_2, color = NMF),
        alpha = 0.8, size = 2
    ) +
    scale_color_discrete() +
    labs(title = "Normal adult mouse projection") +
    theme(legend.position = "bottom")
ggsave(here("output/integ_proj_nmf.png"), plot = combined_plot, width = 10, height = 5)

# Heatmap
plot_data <- as.data.frame.matrix(table(tumor$predicted.id, tumor$NMF)) |>
    tibble::rownames_to_column("predicted_id") |>
    tidyr::pivot_longer(-predicted_id,
        names_to = "NMF",
        values_to = "count"
    )

hm <- ggplot(plot_data, aes(x = NMF, y = predicted_id, fill = count)) +
    geom_tile() +
    geom_text(
        aes(label = ifelse(count > 0, scales::comma(count), "")),
        color = "white",
        size = 3
    ) +
    scale_fill_viridis_c(
        option = "plasma",
        direction = -1,
        name = "Count"
    ) +
    theme_minimal() +
    labs(
        title = "Heatmap of Predicted IDs vs. NMF Modules",
        x = "NMF Module",
        y = "Predicted ID",
        fill = "Count"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
    )

ggsave(here("output/integ_heatmap.png"), plot = hm, width = 10, height = 5)

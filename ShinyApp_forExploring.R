library(shiny)
library(Seurat)
library(ggplot2)
library(patchwork)

# Increase file upload size limit
options(shiny.maxRequestSize = 200*1024^2)  # Set to 200 MB

# Define UI
ui <- fluidPage(
  titlePanel("10X Genomics Data Integration and Visualization"),
  sidebarLayout(
    sidebarPanel(
      fileInput("controlData", "Upload Control Data (filtered_feature_bc_matrix.h5)", accept = c(".h5")),
      fileInput("treatedData", "Upload Treated Data (filtered_feature_bc_matrix.h5)", accept = c(".h5")),
      actionButton("process", "Run Analysis"),
      selectInput(
        "plotType",
        "Select Visualization Type",
        choices = c("QC Violin Plots", "UMAP: Before Integration", "UMAP: After Integration", "Cluster Distribution"),
        selected = "QC Violin Plots"
      )
    ),
    mainPanel(
      plotOutput("mainPlot"),
      tableOutput("clusterCounts")
    )
  )
)

# Define Server
server <- function(input, output, session) {
  observeEvent(input$process, {
    req(input$controlData, input$treatedData)
    
    # Load data
    control_path <- input$controlData$datapath
    treated_path <- input$treatedData$datapath
    
    Control.data <- Read10X_h5(filename = control_path)
    NVP2.data <- Read10X_h5(filename = treated_path)
    
    # Create Seurat objects
    Control.data <- CreateSeuratObject(counts = Control.data, project = "Control")
    NVP2.data <- CreateSeuratObject(counts = NVP2.data, project = "Treated")
    
    # QC metrics
    NVP2.data[["percent.mt"]]  <- PercentageFeatureSet(NVP2.data, pattern = "^MT-")
    NVP2.data[["percent.rbp"]] <- PercentageFeatureSet(NVP2.data, pattern = "^RP[SL]")
    Control.data[["percent.mt"]]  <- PercentageFeatureSet(Control.data, pattern = "^MT-")
    Control.data[["percent.rbp"]] <- PercentageFeatureSet(Control.data, pattern = "^RP[SL]")
    
    # Filtering
    NVP2.data <- subset(NVP2.data, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
    Control.data <- subset(Control.data, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
    
    # Integration preparation
    CN_list <- list(NVP2.data = NVP2.data, Control.data = Control.data)
    for (i in seq_along(CN_list)) {
      CN_list[[i]] <- NormalizeData(CN_list[[i]], verbose = FALSE)
      CN_list[[i]] <- FindVariableFeatures(CN_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    }
    
    # Integration
    CN_anchors <- FindIntegrationAnchors(object.list = CN_list, dims = 1:30)
    CN_seurat <- IntegrateData(anchorset = CN_anchors, dims = 1:30)
    
    # Un-integrated workflow
    DefaultAssay(CN_seurat) <- "RNA"
    CN_seurat <- NormalizeData(CN_seurat, verbose = TRUE)
    CN_seurat <- FindVariableFeatures(CN_seurat, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
    CN_seurat <- ScaleData(CN_seurat, verbose = TRUE)
    CN_seurat <- RunPCA(CN_seurat, npcs = 30, verbose = TRUE)
    CN_seurat <- RunUMAP(CN_seurat, reduction = "pca", dims = 1:30, verbose = TRUE)
    
    # Integrated workflow
    DefaultAssay(CN_seurat) <- "integrated"
    CN_seurat <- ScaleData(CN_seurat, verbose = TRUE)
    CN_seurat <- RunPCA(CN_seurat, npcs = 30, verbose = TRUE)
    CN_seurat <- RunUMAP(CN_seurat, reduction = "pca", dims = 1:30, verbose = TRUE)
    
    # Clustering
    CN_seurat <- FindNeighbors(CN_seurat, dims = 1:30, k.param = 10, verbose = FALSE)
    CN_seurat <- FindClusters(CN_seurat, verbose = FALSE)
    cluster_counts <- table(CN_seurat@meta.data$seurat_clusters, CN_seurat@meta.data$orig.ident)
    
    # Cell cycle heatmap
    cc.genes <- cc.genes.updated.2019
    cell_cycle_genes <- unique(c(cc.genes$s.genes, cc.genes$g2m.genes))
    cell_cycle_genes <- intersect(cell_cycle_genes, rownames(CN_seurat))
    CN_seurat <- ScaleData(CN_seurat, features = cell_cycle_genes)
    
    # Render plots
    output$mainPlot <- renderPlot({
      if (input$plotType == "QC Violin Plots") {
        treated_plot <- VlnPlot(NVP2.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), ncol = 4) +
          ggtitle("Treated Sample QC Metrics")
        control_plot <- VlnPlot(Control.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), ncol = 4) +
          ggtitle("Control Sample QC Metrics")
        treated_plot | control_plot
      } else if (input$plotType == "UMAP: Before Integration") {
        DimPlot(CN_seurat, reduction = "umap", group.by = "orig.ident") + 
          ggtitle("Control and Treated Cells: Before Integration")
      } else if (input$plotType == "UMAP: After Integration") {
        DimPlot(CN_seurat, reduction = "umap", group.by = "orig.ident") + 
          ggtitle("Control and Treated Cells: After Integration")
      } else if (input$plotType == "Cluster Distribution") {
        DimPlot(CN_seurat, reduction = "umap", label = TRUE) +
          ggtitle("Clustered Cells")
      } else if (input$plotType == "Cell Cycle Heatmap") {
        DoHeatmap(CN_seurat, features = cell_cycle_genes, group.by = "orig.ident") +
          ggtitle("Cell Cycle Genes Expression")
      }
    })
    
    # Render cluster count table
    output$clusterCounts <- renderTable({
      as.data.frame(cluster_counts)
    })
  })
}

# Run the app
shinyApp(ui = ui, server = server)

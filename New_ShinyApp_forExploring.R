library(shiny)
library(Seurat)
library(ggplot2)
library(patchwork)


# Increase file upload size limit
options(shiny.maxRequestSize = 200*1024^2)  # Set to 200 MB

# Define UI
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "darkly"),
  titlePanel("10X Genomics Data Integration and Visualization"),
  tabsetPanel(
    tabPanel("Import Data", 
             fileInput("controlData", "Upload Control Data (filtered_feature_bc_matrix.h5)", accept = c(".h5")),
             fileInput("treatedData", "Upload Treated Data (filtered_feature_bc_matrix.h5)", accept = c(".h5")),
             actionButton("process", "Run QC & Integration"),
             selectInput("plotType", "Select QC Plot Type:", 
                         choices = c("QC Violin Plots", "Feature Scatter"), selected = "QC Violin Plots"),
             plotOutput("qcPlot"),
             downloadButton("downloadQCPlot", "Download QC Plot")
    ),
    tabPanel("Integration Analysis", 
             plotOutput("umapPlot"),
             downloadButton("downloadUMAP", "Download UMAP Plot"),
             plotOutput("clusterPlot"),
             downloadButton("downloadClusterPlot", "Download Cluster Plot"),
             tableOutput("clusterTable"),
             downloadButton("downloadClusterTable", "Download Cluster Data")
    ),
    tabPanel("Cell Cycle Analysis", 
             plotOutput("cellCyclePlot"),
             downloadButton("downloadCellCyclePlot", "Download Cell Cycle Heatmap")
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  # Reactive values to store Seurat objects
  seurat_objects <- reactiveValues(control = NULL, treated = NULL, integrated = NULL)
  
  observeEvent(input$process, {
    req(input$controlData, input$treatedData)
    
    # Load data
    control_path <- input$controlData$datapath
    treated_path <- input$treatedData$datapath
    
    Control.data <- Read10X_h5(filename = control_path)
    Treated.data <- Read10X_h5(filename = treated_path)
    
    # Create Seurat objects
    Control.data <- CreateSeuratObject(counts = Control.data, project = "Control")
    Treated.data <- CreateSeuratObject(counts = Treated.data, project = "Treated")
    
    # Compute QC metrics
    Control.data[["percent.mt"]]  <- PercentageFeatureSet(Control.data, pattern = "^MT-")
    Control.data[["percent.rbp"]] <- PercentageFeatureSet(Control.data, pattern = "^RP[SL]")
    Treated.data[["percent.mt"]]  <- PercentageFeatureSet(Treated.data, pattern = "^MT-")
    Treated.data[["percent.rbp"]] <- PercentageFeatureSet(Treated.data, pattern = "^RP[SL]")
    
    # **Filtering**
    Treated.data <- subset(Treated.data, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
    Control.data <- subset(Control.data, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)
    
    # **Integration Preparation**
    CN_list <- list(Treated.data = Treated.data, Control.data = Control.data)
    for (i in seq_along(CN_list)) {
      CN_list[[i]] <- NormalizeData(CN_list[[i]], verbose = FALSE)
      CN_list[[i]] <- FindVariableFeatures(CN_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    }
    
    # **Integration**
    CN_anchors <- FindIntegrationAnchors(object.list = CN_list, dims = 1:30)
    CN_seurat <- IntegrateData(anchorset = CN_anchors, dims = 1:30)
    
    # **Un-integrated Workflow**
    DefaultAssay(CN_seurat) <- "RNA"
    CN_seurat <- NormalizeData(CN_seurat, verbose = TRUE)
    CN_seurat <- FindVariableFeatures(CN_seurat, selection.method = "vst", nfeatures = 2000, verbose = TRUE)
    CN_seurat <- ScaleData(CN_seurat, verbose = TRUE)
    CN_seurat <- RunPCA(CN_seurat, npcs = 30, verbose = TRUE)
    CN_seurat <- RunUMAP(CN_seurat, reduction = "pca", dims = 1:30, verbose = TRUE)
    
    # **Integrated Workflow**
    DefaultAssay(CN_seurat) <- "integrated"
    CN_seurat <- ScaleData(CN_seurat, verbose = TRUE)
    CN_seurat <- RunPCA(CN_seurat, npcs = 30, verbose = TRUE)
    CN_seurat <- RunUMAP(CN_seurat, reduction = "pca", dims = 1:30, verbose = TRUE)
    
    # **Clustering**
    CN_seurat <- FindNeighbors(CN_seurat, dims = 1:30, k.param = 10, verbose = FALSE)
    CN_seurat <- FindClusters(CN_seurat, verbose = FALSE)
    
    # Store integrated object
    seurat_objects$control <- Control.data
    seurat_objects$treated <- Treated.data
    seurat_objects$integrated <- CN_seurat
  })
  
  # **QC Plots**
  output$qcPlot <- renderPlot({
    req(seurat_objects$control, seurat_objects$treated)
    
    if (input$plotType == "QC Violin Plots") {
      treated_plot <- VlnPlot(seurat_objects$treated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), ncol = 4) +
        ggtitle("Treated Sample QC Metrics")
      control_plot <- VlnPlot(seurat_objects$control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"), ncol = 4) +
        ggtitle("Control Sample QC Metrics")
      treated_plot | control_plot
    }else if (input$plotType == "Feature Scatter") {
      treated_plot <- FeatureScatter(seurat_objects$treated, feature1 = "nCount_RNA", feature2 = "percent.mt") +
        ggtitle("Treated Sample Feature Scatter")
      control_plot <- FeatureScatter(seurat_objects$control, feature1 = "nCount_RNA", feature2 = "percent.mt") +
        ggtitle("Control Sample Feature Scatter")
      treated_plot | control_plot
    }
    
  })
  
  output$downloadQCPlot <- downloadHandler(
    filename = "QC_Plot.png",
    content = function(file) {
      ggsave(file, plot = output$qcPlot())
    }
  )
  
  # **UMAP Plot**
  output$umapPlot <- renderPlot({
    req(seurat_objects$integrated)
    DimPlot(seurat_objects$integrated, reduction = "umap", group.by = "seurat_clusters") +
      ggtitle("UMAP of Integrated Data")
  })
  
  output$downloadUMAP <- downloadHandler(
    filename = "UMAP_Plot.png",
    content = function(file) {
      ggsave(file, plot = output$umapPlot())
    }
  )
  
  # **Cluster Plot**
  output$clusterPlot <- renderPlot({
    req(seurat_objects$integrated)
    FeaturePlot(seurat_objects$integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), reduction = "umap") +
      ggtitle("Feature Distribution Across Clusters")
  })
  
  output$downloadClusterPlot <- downloadHandler(
    filename = "Cluster_Plot.png",
    content = function(file) {
      ggsave(file, plot = output$clusterPlot())
    }
  )
  
  # **Cluster Table**
  output$clusterTable <- renderTable({
    req(seurat_objects$integrated)
    table(seurat_objects$integrated@meta.data$seurat_clusters, seurat_objects$integrated@meta.data$orig.ident)
  })
  
  output$downloadClusterTable <- downloadHandler(
    filename = "Cluster_Assignments.csv",
    content = function(file) {
      write.csv(output$clusterTable(), file)
    }
  )
  
  # **Cell Cycle Plot**
  output$cellCyclePlot <- renderPlot({
    req(seurat_objects$integrated)
    DoHeatmap(seurat_objects$integrated, features = cc.genes.updated.2019$s.genes, group.by = "seurat_clusters") +
      ggtitle("Cell Cycle Gene Expression")
  })
  
  output$downloadCellCyclePlot <- downloadHandler(
    filename = "Cell_Cycle_Heatmap.png",
    content = function(file) {
      ggsave(file, plot = output$cellCyclePlot())
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)

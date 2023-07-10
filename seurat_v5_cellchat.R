library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(Azimuth)
library(BPCells)
options(Seurat.object.assay.version = "v5")
parse.data <- open_matrix_anndata_hdf5(
"/home/baranov_lab/10X/HumanCellAtlas/raw/snRNA_rawcounts.h5ad"
)
write_matrix_dir(mat = parse.data, dir = "/home/baranov_lab/10X/HumanCellAtlas/raw/parse_1m_pbmc")
parse.mat <- open_matrix_dir(dir = "/home/baranov_lab/10X/HumanCellAtlas/raw/parse_1m_pbmc")
library(data.table)
metadata <- fread("/home/baranov_lab/10X/HumanCellAtlas/raw-metadata.csv",header=TRUE)
head(metadata)
head(parse.mat)
rownames(metadata) <- metadata$V1
head(metadata)
TOTAL <- CreateSeuratObject(counts = parse.mat, meta.data = metadata)
head(TOTAL)
TOTAL <- NormalizeData(TOTAL)
TOTAL[['RNA']] <- split(TOTAL[['RNA']], f = TOTAL$majorclass)
TOTAL <- FindVariableFeatures(TOTAL, nfeatures = 3000)
TOTAL <- SketchData(TOTAL, ncells = 50000, method = 'LeverageScore', sketched.assay = 'sketch_split') ##SHOULD BE RAW DATA
DefaultAssay(TOTAL) <- 'sketch_split'
TOTAL <- FindVariableFeatures(TOTAL, nfeatures = 3000)
TOTAL <- ScaleData(TOTAL)
TOTAL <- RunPCA(TOTAL)
TOTAL <- IntegrateLayers(TOTAL, method = RPCAIntegration, orig = 'pca', new.reduction = 'integrated.rpca', dims = 1:30, k.anchor = 20, 
                         reference = which(Layers(TOTAL, search = 'data') %in% c('data.Rod')))
TOTAL <- FindNeighbors(TOTAL, reduction = 'integrated.rpca', dims = 1:30)
TOTAL <- FindClusters(TOTAL, resolution = 1)
TOTAL <- RunUMAP(TOTAL, reduction = 'integrated.rpca', dims = 1:30, return.model = T)
TOTAL[['sketch_split']] <- JoinLayers(TOTAL[['sketch_split']])
total_sketch <- DietSeurat(TOTAL, assays = 'sketch_split')
saveRDS(total_sketch, '/home/baranov_lab/10X/HumanCellAtlas/raw/TOTAL_sketch_split300k.rds')
gc()

library(CellChat)

cellchat <- createCellChat(object = total_sketch, group.by = "majorclass")
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, trim = 0.25, population.size = FALSE, raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 50)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 30)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height  = 30)
ht1 + ht2
saveRDS(cellchat, '/home/baranov_lab/10X/HumanCellAtlas/raw/cellchat_sketch_split300k.rds')


library(ArchR)
library(parallel)
set.seed(42)

# Get args
args = commandArgs(trailingOnly=TRUE)
sample_name = args[1]
cr_outs = args[2]
out_dir = args[3]

# Configure
addArchRThreads(threads = 25) 
addArchRGenome('hg38')

# Make a new ArchR project for each sample
current_dir <- getwd()
print(paste0("Processing ", sample_name))

# set input path
input_files <- c(file.path(cr_outs,'atac_fragments.tsv.gz'))
names(input_files) <- c(sample_name)

# set output directory
output_dir <- file.path(out_dir,sample_name, 'ArchR')
dir.create(output_dir)
setwd(output_dir)

#create arrow files
print(paste0("Creating ArrowFiles "))
ArrowFiles <- createArrowFiles(
    inputFiles = c(input_files),
    sampleNames = names(input_files),
    minTSS = 0,
    minFrags = 1, 
    addTileMat = TRUE,
    maxFrags = 1e+20,
    addGeneScoreMat = FALSE,
    excludeChr = c('chrM'))
    
# create project
addArchRThreads(threads = 1) 
proj_name <- 'pre_filt'
print(paste0("Creating ArchR project. Project name:",proj_name))
proj <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = proj_name,
    copyArrows = FALSE
)
    
# save project
print(paste0("Saving ArchR project. Project name:",proj_name))
proj <- saveArchRProject(ArchRProj = proj)
    
# export metadata
dir.create(file.path(proj_name, 'export'))
print(paste0("Exporting cell_metadata.csv"))
write.csv(getCellColData(proj), 
            file.path(proj_name, 'export', 'cell_metadata.csv'), 
            quote=FALSE)
    
# get cellranger filtered cells
cr_barcodes <- read.table(file.path(cr_outs,'filtered_feature_bc_matrix/barcodes.tsv.gz'),
                            header = FALSE,
                            col.names=c('barcode'))
cr_barcodes$barcode <- paste0(sample_name,'#',cr_barcodes$barcode)
    
# subset project
proj_name <- 'cr_filt'
print(paste0("Creating ArchR project. Project name:",proj_name))
sub_proj <- subsetArchRProject(proj,
                                cells = cr_barcodes$barcode,
                                outputDirectory = proj_name,
                                dropCells = TRUE)
    
# run iterative LSI and UMAP
print("Running Itterative LSI and adding UMAP")
sub_proj <- addIterativeLSI(ArchRProj = sub_proj, useMatrix = "TileMatrix", name = "IterativeLSI")
sub_proj <- addUMAP(sub_proj)
sub_proj <- saveArchRProject(ArchRProj = sub_proj)
print(paste0("Saving ArchR project. Project name:",proj_name))
    
# export metadata
print("Exporting cell_metadata.csv")
write.csv(getCellColData(sub_proj), 
            file.path(proj_name, 'export', 'cell_metadata.csv'), 
            quote=FALSE)
print("Exporting svd.csv")
write.csv(getReducedDims(sub_proj), 
            file.path(proj_name, 'export', 'svd.csv'), 
            quote=FALSE)
print("Exporting umap.csv")
write.csv(getEmbedding(sub_proj), 
            file.path(proj_name, 'export', 'umap.csv'), 
            quote=FALSE)
    
#reset working directory
setwd(current_dir)

library(ArchR)
library(parallel)
set.seed(42)

# Get args
args = commandArgs(trailingOnly=TRUE)
data_dir = args[1]
out_dir = args[2]

# Configure
addArchRThreads(threads = 25) 
addArchRGenome('hg38')

# Make a new ArchR project for each sample
sample_names <- list.dirs(data_dir,recursive = FALSE, full.names = FALSE)
current_dir <- getwd()
for (sample_name in sample_names) {
    print(paste0("Processing ", sample_name))
    # set input path
    input_files <- c(file.path(data_dir,sample_name,'atac_fragments.tsv.gz'))
    names(input_files) <- c(sample_name)

    # set output directory
    output_dir <- file.path(out_dir,sample_name, 'ArchR')
    dir.create(output_dir)
    setwd(output_dir)

    #create arrow files
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
    proj <- ArchRProject(
        ArrowFiles = ArrowFiles, 
        outputDirectory = sample_name,
        copyArrows = FALSE
    )
    
    # save project
    proj <- saveArchRProject(ArchRProj = proj)
    
    # export metadata
    dir.create(file.path(sample_name, 'export'))
    print(file.path(sample_name, 'export'))
    write.csv(getCellColData(proj), 
              file.path(sample_name, 'export', 'cell_metadata.csv'), 
              quote=FALSE)
    
    #reset working directory
    setwd(current_dir)
    
    }

gene_list <- function(file, path = "/Users/nacho/gene_info/gene_sets"){
    genes <- scan(file.path(path, file), "character")
    genes
}

genes_from_coords <- gfc <- function(coords){
    sapply(coords, function(x) str_match(x,"(^.+)_[^_]+_[^_]+$")[2]) # we have to account for "A_33_P3254580_28_29"
}

coords_from_genes <- cfg <- function(gene_names, genes){
    rows <- genes[genes$GeneName %in% gene_names, ]
    rownames(rows)
}

validate_paths <- function(paths) {
    sapply(paths, function(path){
        if(!file.exists(path)) {
            dir.create(path, recursive=T)
            message("Created path ", path)
        }
    })
    invisible()
}

fix_gene_symbols <- function(genes, translator_path){
    agilent_translator <- read.table(translator_path, sep=",", stringsAsFactors=F)
    colnames(agilent_translator) <- c("ProbeName","GeneName")

    # extract the ProbeName and GeneName of the genes that have translations
    from <- unique(subset(genes[genes$ProbeName %in% agilent_translator$ProbeName,], select=c(ProbeName,GeneName)))
    from <- from[order(from$ProbeName),]

    # make sure translator matches the order and the probe names in "from" (we call this subset "to")
    to <- data.frame(ProbeName=from$ProbeName, GeneName=agilent_translator$GeneName[match(from$ProbeName, agilent_translator$ProbeName)], stringsAsFactors=FALSE)

    if (nrow(to) != nrow(from)){
        stop("Repeated probe names in agilent translator file ", translator_path, "'from' and 'to' have different number of rows")
    }

    translatable_positions <- which(genes$GeneName %in% from$GeneName)
    translated_gene_names <- to$GeneName[match(genes$GeneName[translatable_positions],from$GeneName)]
    genes$GeneName[translatable_positions] <- translated_gene_names

    genes
}


get_RG <- function(d){
    message("Generating RG")
    message("Cy3: ", paste(unique(d$targets$Cy3), collapse = " "))
    message("Cy5: ", paste(unique(d$targets$Cy5), collapse = " "))

    RG <- read.maimages(d$targets, source = "agilent", path = d$paths$raw_data)
    RG$old_genes <- RG$genes

    message("Updating obsolete gene symbols")
        RG$genes <- fix_gene_symbols(RG$genes, d$paths$agilent_translator)
        rownames(RG$genes) <- paste0(RG$genes$GeneName, "_", RG$genes$Row, "_", RG$genes$Col)
    message(length(RG$genes$GeneName[RG$genes$GeneName != RG$old_genes$GeneName])," gene symbols fixed") # this is an overestimation because some might be duplicates

    for (obj in list("G", "Gb", "R", "Rb")){
        colnames(RG[[obj]]) <- d$targets$Id
        rownames(RG[[obj]]) <- rownames(RG$genes)
    }

    saveRDS(RG, file = file.path(d$paths$calculations, paste0(d$name, "_RG.rds")))

    RG
}

get_MAwin <- function(d){
    message("Generating MAwin")

    MAwin <- normalizeWithinArrays(d$RG, method="loess", bc.method="subtract")
    rownames(MAwin$M) <- paste0(MAwin$genes$GeneName, "_", MAwin$genes$Row, "_", MAwin$genes$Col)

    # remove control probes
    MAwin <- MAwin[MAwin$genes$ControlType == 0,]

    saveRDS(MAwin, file = file.path(d$paths$calculations, paste0(d$name, "_MAwin.rds")))

    MAwin
}

# column_indeces is a vector that assigns a column to be substracted to each column.
# ex: m1_d0 m1_d3 m1_d5 m2_d0 m2_d4 => 1 1 1 4 4 (columns 1 to 3 - column 1, columns 4 to 5 - column 4)
get_M_zero <- function(d, column_indeces){
    M_zero <- d$M - d$M[, column_indeces]
    M_zero
}

get_design_by_serotype <- function(targets) {

    initial_challenge <- factor(targets$initial_challenge)
    time <- factor(targets$Cy3)
    initial_challenge_time <- factor(paste(initial_challenge, time, sep = "_"), levels = unique(paste(initial_challenge, time, sep = "_")))
    design <- model.matrix(~ 0 + initial_challenge_time)

    colnames(design) <- as.character(levels(initial_challenge_time))

    design
}

get_fit <- function(d, m_values = d$M, design = d$design){
    message("Generating fit")

    fit <- lmFit(m_values, design)
    fit <- eBayes(fit)

    fit
}

get_fit_correlation <- function(d, m_values = d$M, design = d$design, block = d$targets$sample_name){
    message("Generating fit")

    corrfit <- duplicateCorrelation(m_values, design, block = block)
    fit <- lmFit(m_values, design, block = block, correlation = corrfit$consensus)
    fit <- eBayes(fit)

    fit
}


get_de <- function(d, contrast, fit = d$fit) {

    cont.matrix <- makeContrasts(contrasts = contrast, levels = fit$coefficients)
    fit <- contrasts.fit(fit, cont.matrix)
    fit <- eBayes(fit) # you need to calculate the eBayes statistic after every contrast

    de <- topTable(fit, coef = 1, adjust.method = "BH", number = nrow(fit))
    de <- de[with(de, order(adj.P.Val)), ] # sort de by ascending p-value
    rownames(de) <- de$ID
    de <- de[,-1]
    de$GeneName <- genes_from_coords(rownames(de))
    de$visual_pvalue <- -log10(de$adj.P.Val)

    de
}





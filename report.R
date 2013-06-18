library(df2json)
library(yaml)
library(rjson)
library(reshape2)
library(plyr)
library(clickme)
library(stringr)
library(statmod)
library(limma)

source("main.r")

d <- list()

d$dengue$name <- "dengue"
d$dengue$info <- yaml.load_file("dengue_targets.yaml")
d$dengue$paths <- d$dengue$info$paths
validate_paths(d$dengue$paths)
d$dengue$targets <- ldply(d$dengue$info$targets, data.frame)
d$dengue$targets$infections <- toupper(paste(d$dengue$targets$initial_challenge, "denv5", sep = " + "))
d$dengue$arrays <- function(timepoint){
    which(d$dengue$targets$Cy3 == timepoint)
}

d$dengue$RG <- get_RG(d$dengue)
d$dengue$genes <- d$dengue$RG$genes

d$dengue$MAwin <- get_MAwin(d$dengue)

d$dengue$M <- d$dengue$MAwin$M

d$dengue$M_zero <- get_M_zero(d$dengue, rep(d$dengue$array(0), each = 4))

d$dengue$mds <- plotMDS(d$dengue$M, gene.selection = "common", top = 2000, labels = d$dengue$targets$group, col = color(d$dengue$targets$group), cex = 2, main = "dengue arrays\nby id")

d$dengue$num_nas <- apply(d$dengue$M, 2, function(x) sum(is.na(x)))

saveRDS(d$dengue, "/Users/nacho/Documents/BU/Connor/projects/dengue/calculations/dengue.rds")

# after the above runs once I usually comment it out and load the rds object directly
d$dengue <- readRDS("/Users/nacho/Documents/BU/Connor/projects/dengue/calculations/dengue.rds")
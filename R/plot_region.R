#!/usr/bin/env Rscript
library(optparse)
library(tidyverse)

options(stringsAsFactors=FALSE)

source("helpers.R")

option_list = list(
    make_option(
        c("-i", "--infiles"), type="character",
        help="Comma-separated list of bedgraph files containing scores to be plotted",
    ),
    make_option(
        c("-s", "--samplenames"), type="character",
        help="Comma-separated list of sample names associated with each of the bedgraph files passed to `--infiles`",
    ),
    make_option(
        c("-f", "--features"), type="character", default=NULL,
        help="Optional. Contains tab-separated feature information for plotting with signals present in `infile`.",
    ),
    make_option(
        c("-c", "--contig"), type="character",
        help="The name of the contig to plot",
    ),
    make_option(
        c("-l", "--start"), type="integer",
        help="The zero-indexed start position of the plot",
    ),
    make_option(
        c("-e", "--end"), type="integer",
        help="The zero-indexed end position of the plot",
    ),
    make_option(
        c("-p", "--plot_file"), type="character",
        help="File created containing the plot",
    ),
    make_option(
        c("--plot_height"), type="double", default=4,
        help="Plot hieght in inches",
    ),
    make_option(
        c("--plot_width"), type="double", default=4,
        help="Plot width in inches",
    ),
    make_option(
        c("--line_alpha"), type="double", default=1.0,
        help="Opacity of lines, 1 is opaque, 0 is transparent",
    ),
    make_option(
        c("--line_size"), type="double", default=1.25,
        help="Scaling factor for lines",
    ),
    make_option(
        c("--yvar"), type="character", default="score",
        help="Name of the y-axis variable",
    ),
    make_option(
        c("--upper_files"), type="character", default=NULL,
        help="Name of the bedgraph files containing the upper CI. Must be in same order as files in --infiles",
    ),
    make_option(
        c("--lower_files"), type="character", default=NULL,
        help="Name of the bedgraph files containing the lower CI. Must be in same order as files in --infiles",
    ),
    make_option(
        c("--colorvar"), type="character", default="sample",
        help="Name of the color variable, if such a variable exists",
    ),
    make_option(
        c("--linetypevar"), type="character", default=NULL,
        help="Name of the linetype variable, if such a variable exists",
    ),
    make_option(
        c("--ylabel"), type="character", default="score",
        help="Label assigned to y-axis",
    ),
    make_option(
        c("--featurevar"), type="character", default="name",
        help="Variable in `features` with feature names",
    ),
    make_option(
        c("--log"), type="logical", action="store_true", default=FALSE,
        help="Include at command line if you want to plot the y-axis as a log scale",
    ),
    make_option(
        c("--facet"), type="character", default=NULL,
        help="Passed to facet_grid internally to set up facetted plotting",
    ),
    make_option(
        c("--name_angle"), type="numeric", default=0,
        help="angle to rotate feature names by."
    ),
    make_option(
        c("--include_feature_types"), type="character", default="CDS,tRNA,rRNA",
        help="Comma-separated list of feature types to include (default is 'CDS,tRNA,rRNA')"
    )
)
 
print("Reading command line options")
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

infiles = str_split(opt$infiles, ",", simplify=TRUE)[1,]
samples = str_split(opt$samplenames, ",", simplify=TRUE)[1,]
feature_types = str_split(opt$include_feature_types, ",", simplify=TRUE)[1,]
up=NULL
low=NULL
upper_files = NULL
lower_files = NULL

if (length(infiles) != length(samples)) {
    stop("The number of input files must match the number of samples")
}

if (!is.null(opt$lower_files)) {
    if (is.null(opt$lower_files)) {
        stop("If you specify lower bounds, you must also specify upper bounds.")
    }
    lower_files = str_split(opt$lower_files, ",", simplify=TRUE)[1,]
    if (length(infiles) != length(lower_files)) {
        stop("The number of input files must match the number of upper limit files")
    }
    up = "upper"
}

if (!is.null(opt$upper_files)) {
    if (is.null(opt$lower_files)) {
        stop("If you specify upper bounds, you must also specify lower bounds.")
    }
    upper_files = str_split(opt$upper_files, ",", simplify=TRUE)[1,]
    if (length(infiles) != length(upper_files)) {
        stop("The number of input files must match the number of upper limit files")
    }
    low = "lower"
}

for (i in 1:length(infiles)) {
    colnames = c("seqname","start","end","score")
    coltypes = c("c","i","i","d")
    if (!is.null(opt$facet)) {
        colnames[5] = str_remove(opt$facet, " ~")
        coltypes[5] = "c"
    }
    print(colnames)
    coltypes = paste(coltypes, collapse="")
    this_file = infiles[i]
    this_sample = samples[i]
    print(paste0("Reading data from ", this_file))
    tmp_tib = read_table(
        file=this_file,
        col_names=colnames,
        col_types=coltypes,
    ) %>% mutate(sample = this_sample)
    if (!is.null(opt$upper_files)) {
        this_upper_file = upper_files[i]
        this_lower_file = lower_files[i]
        tmp_upper = read_table(
            file=this_upper_file,
            col_names=colnames,
            col_types=coltypes,
        ) %>% mutate(sample = this_sample)
        tmp_lower = read_table(
            file=this_lower_file,
            col_names=colnames,
            col_types=coltypes,
        ) %>% mutate(sample = this_sample)
        tmp_tib = tmp_tib %>%
            mutate(upper=tmp_upper$score, lower=tmp_lower$score)
    }

    if (i == 1) {
        data_tib = tmp_tib
    } else {
        data_tib = data_tib %>% bind_rows(tmp_tib)
    }
}

if (!is.null(opt$features)) {
    plot_features = TRUE
    gff = gffRead(
        opt$features
    )  %>% filter(
        seqname == opt$contig,
        feature %in% feature_types
    )
    gff$name = getAttributeField(gff$attributes, "Name")
    gff$locus_tag = getAttributeField(gff$attributes, "locus_tag")
    gff = gff %>% mutate(
        name = ifelse(duplicated(name), locus_tag, name)
    )
} else {
    plot_features = FALSE
    gff = NULL
}

p = plot_locus(
    data_tib,
    opt$start,
    opt$end,
    opt$contig,
    feats_df = gff,
    lineAlpha = opt$line_alpha,
    lineSize = opt$line_size,
    yvar = opt$yvar,
    color_var = opt$colorvar,
    linetype_var = opt$linetypevar,
    ylabel = opt$ylabel,
    feat_var = opt$featurevar,
    facet = opt$facet,
    plotFeatures = plot_features,
    log_y = opt$log,
    feat_fill_var = "feature",
    name_angle = opt$name_angle,
    upper = up,
    lower = low
)
ggsave(plot=p, filename=opt$plot_file, height=opt$plot_height, width=opt$plot_width)


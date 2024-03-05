#!/usr/bin/env Rscript
library(optparse)
library(tidyverse)
library(scales)

get_exec_file = function() {
    # copied from https://stackoverflow.com/questions/1815606/determine-path-of-the-executing-script
    cmdArgs = commandArgs(trailingOnly = FALSE)
    needle = "--file="
    match = grep(needle, cmdArgs)
    if (length(match) > 0) {
        # Rscript
        return(normalizePath(sub(needle, "", cmdArgs[match])))
    } else {
        # 'source'd via R console
        return(normalizePath(sys.frames()[[1]]$ofile))
    }
}

options(stringsAsFactors=FALSE)

this_file = get_exec_file()
this_direc = dirname(this_file)
source(file.path(this_direc,"helpers.R"))

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
        help="Optional. A gff file with feature annotations.",
    ),
    make_option(
        c("--feature_name_field"), type="character", default=NULL,
        help="Optional. Sets which field in the attribute column of your gff file contains the feature names to be plotted. If not provided, will use, in order of preference, gene, Name, locus_tag, and ID any of those fields exist. Note that for features to plot correctly, the chosen field MUST be uniqe to each feature to be plotted.",
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
        c("--colorvalues"), type="character", default=NULL,
        help="Optional comma-separated list of hexadecimal colors. If provided user must set a color value for each DISTINCT value provided in --samplenames.",
    ),
    make_option(
        c("--linetypevar"), type="character", default=NULL,
        help="Name of the linetype variable, if such a variable exists"
    ),
    make_option(
        c("--ylabel"), type="character", default="score",
        help="Label assigned to y-axis"
    ),
    make_option(
        c("--ylims"), type="character", default=NULL,
        help="Limits (using coord_cartesian) to place on y-axis"
    ),
    make_option(
        c("--ygrid"), type="logical", action="store_true", default=FALSE,
        help="Include at command line if you want to plot a grid of horizonal lines."
    ),
    make_option(
        c("--featurevar"), type="character", default="name",
        help="Variable in `features` with feature names"
    ),
    make_option(
        c("--log"), type="logical", action="store_true", default=FALSE,
        help="Include at command line if you want to plot the y-axis as a log scale"
    ),
    make_option(
        c("--facet"), type="character", default=NULL,
        help="Passed to facet_grid internally to set up facetted plotting"
    ),
    make_option(
        c("--name_angle"), type="numeric", default=0,
        help="angle to rotate feature names by."
    ),
    make_option(
        c("--include_feature_types"), type="character", default="CDS,tRNA,rRNA",
        help="Comma-separated list of feature types to include (default is 'CDS,tRNA,rRNA')"
    ),
    make_option(
        c("--feature_colors"), type="character", default=NULL,
        help="Optional comma-separated list of hexadecimal values assigning a color to each feature type provided by the --include_feature_types argument. If ommitted, ggplot defaults will be used."
    ),
    make_option(
        c("--motifs_file"), type="character", default=NULL,
        help="Optional gff file containing annotations that will be plotted as a lightly-shaded rectangle underneath the data."
    )
)
 
print("Reading command line options")
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

infiles = str_split(opt$infiles, ",", simplify=TRUE)[1,]
samples = str_split(opt$samplenames, ",", simplify=TRUE)[1,]
distinct_samples = unique(samples)
feature_types = str_split(opt$include_feature_types, ",", simplify=TRUE)[1,]
if (opt$include_feature_types == "CDS,tRNA,rRNA") {
    if (is.null(opt$feature_colors)) {
        feature_colors = c("#43BA37", "#619CFF", "#F0766D")
    } else {
        feature_colors = str_split(opt$feature_colors, ",", simplify=TRUE)[1,]
    }
} else {
    if (is.null(opt$feature_colors)) {
        # get ggplot default colors
        feature_colors = hue_pal()(length(feature_types))
    } else {
        feature_colors = str_split(opt$feature_colors, ",", simplify=TRUE)[1,]
    }
}

print("ylab argument")
print(opt$ylabel)
ylims = opt$ylims
if (!(is.null(ylims))) {
    ylims = eval(parse(text=ylims))
}

if (length(feature_colors) != length(feature_types)) {
    stop("The number of features provided by --feature_types and the number of hexadecimal colors provided by --feature_colors does not match. Edit the command so they are the same length. Exiting now.")
}

color_vals = NULL
up=NULL
low=NULL
upper_files = NULL
lower_files = NULL

if (length(infiles) != length(samples)) {
    stop("The number of input files must match the number of samples")
}

if (!is.null(opt$colorvalues)) {
    color_vals = str_split(opt$colorvalues, ",", simplify=TRUE)[1,]
    if (length(distinct_samples) != length(color_vals)) {
        print("The number of distinct sample names provided by --samplenames argument must match the number of hexadecimal color values supplied by --colorvalues, but the numbers do not match. For instance, if you provided --samplesnames A,A,B,B,C,C then there are three distinct samples, A,B,C, so you should have three values provided to --colorvalues. ")
        print(paste0("The color_vals you provided were: ", paste(color_vals, collapse=", ")))
        print(paste0("The distinct samples you provided were: ", paste(distinct_samples, collapse=", ")))
        stop("Exiting now.")
    }
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

if (!is.null(opt$motifs_file)) {
    plot_motifs = TRUE
    motif_df = gffRead(
        opt$motifs_file
    ) %>% filter(
        seqname == opt$contig
    )
} else {
    plot_motifs = FALSE
    motif_df = NULL
}

print("================ features arg ===============")
print(opt$features)

if (!is.null(opt$features)) {
    plot_features = TRUE
    gff = gffRead(
        opt$features
    )  %>% filter(
        seqname == opt$contig,
        feature %in% feature_types
    )
    print("========= read gff ============")
    print(head(gff))
    print("========= feature_name_field ============")
    print(opt$feature_name_field)
    if (is.null(opt$feature_name_field)) {
        gff$name = getAttributeField(gff$attributes, "gene")
        gff$alt_name = getAttributeField(gff$attributes, "Name")
        gff$locus_tag = getAttributeField(gff$attributes, "locus_tag")
        gff$ID = getAttributeField(gff$attributes, "ID")
        gff = gff %>% mutate(
            name = ifelse(duplicated(name), alt_name, name)
        )
        gff = gff %>% mutate(
            name = ifelse(duplicated(name), locus_tag, name)
        )
        gff = gff %>% mutate(
            name = ifelse(duplicated(name), ID, name)
        )
        fields = c("gene", "Name", "locus_tag", "ID")
    } else {
        gff$name = getAttributeField(gff$attributes, opt$feature_name_field)
        fields = c(opt$feature_name_field)
    }
} else {
    plot_features = FALSE
    gff = NULL
}

print("========= named gff ============")
print(head(gff))

duped = duplicated(gff$name)

if (any(duped)) {
    dup_names = gff$name[duped]
    dup_idxs = which(duped)
    for (i in 1:length(dup_idxs)) {
        gff$name[dup_idxs[i]] = paste0(gff$name[dup_idxs[i]], "_", i)
    }
    #stop(paste0("Duplicated feature names exist in your gff file! You chose ", fields, " as your gff attribute field(s) to use for feature names.\nDuplicated names:\n", dup_names, "\nEither rename the features so the field is unique to each feature, or consider using a different field altogether."))
}

print("========= dedup gff ============")
print(head(gff))

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
    color_vals = color_vals,
    linetype_var = opt$linetypevar,
    ylabel = opt$ylabel,
    ylims = ylims,
    feat_var = opt$featurevar,
    facet = opt$facet,
    plotFeatures = plot_features,
    feat_types = feature_types,
    feat_colors = feature_colors,
    log_y = opt$log,
    feat_fill_var = "feature",
    name_angle = opt$name_angle,
    ygrid = opt$ygrid,
    upper = up,
    lower = low,
    plotMotifLocs = plot_motifs,
    motifs_df = motifs_df
)

print(paste("Writing plot to ", opt$plot_file))
ggsave(plot=p, filename=opt$plot_file, height=opt$plot_height, width=opt$plot_width)


#!/usr/bin/env Rscript
library(tidyverse)
library(ggnewscale)

options(stringsAsFactors=FALSE)

get_min = function(.df,.contig) {
    return(min(.df$coverage[.df$contig == .contig]))
}

get_max = function(.df, .contig) {
    return(max(.df$coverage[.df$contig == .contig]))
}

get_cor = function(.df, .contig1, .contig2) {
    x = .df$coverage[.df$contig == .contig1]
    y = .df$coverage[.df$contig == .contig2]
    return (cor(x,y,method="pearson"))
}

get_slope = function(.df, .contig1, .contig2) {
    form = paste0(.contig2, " ~ ", .contig1)
    fit = lm(form, .df)
    print(fit$coefficients)
    return (fit$coefficients[2])
}

make_plot = function(.df, .xvar, .yvar, .colvar, .facetvar=NULL) {
    x=sym(.xvar)
    y=sym(.yvar)
    col=sym(.colvar)
    fill=sym(.colvar)
    #facet=sym(.facetvar)
    p = .df %>%
        ggplot(aes(x=!!x, y=!!y, color=!!col, fill=!!col)) +
        #facet_grid(cols=vars(Feature)) +
        geom_point() +
        geom_smooth(method="lm") +
        theme_classic() +
        labs(x=paste0("Spikein: ", .xvar), y="Bacillus chromosome")
    return(p)
}

make_polygons = function(
        row,
        shift_feats,
        shift_ncRNA,
        feat_var,
        plotStart,
        plotEnd,
        shiftTop,
        shiftBottom,
        geneTop,
        geneBottom
) {
    if((shift_ncRNA) && (row$poly_fill %in% shift_feats)) {
        out_tib = create_feature_polygon(
            row$poly_name, row$start, row$end, row$strand,
            top=shiftTop, bottom=shiftBottom,
            plot_width=plotEnd-plotStart
        )
    } else {
        out_tib = create_feature_polygon(
            row$poly_name, row$start, row$end, row$strand,
            top=geneTop, bottom=geneBottom,
            plot_width=plotEnd-plotStart
        )
    }
    return(out_tib)
}

create_feature_polygon = function(name, start, end, strand, top, bottom, plot_width) {
 
    feat_width = end-start
    # start with traingle at 10% of feature width
    #triangle_width = 0.1 * feat_width
    triangle_width = 0.02 * plot_width
    if (triangle_width > 0.5 * feat_width) {
        triangle_width = 0.5 * feat_width
    }
    # make sure triangles cannot be greather than 1-20th the entire width of the plot
    #if (triangle_width > 0.05 * plot_width) {
    #    triangle_width = 0.05 * plot_width
    #}
    # if the triangle will be very narrow, just plot as non-stranded
    if (triangle_width < 0.001 * plot_width) {
        strand = "*"
    }
    if (strand == "+") {
        shape_x = c(start, end-triangle_width, end, end-triangle_width, start)
        shape_y = c(top, top, mean(c(top,bottom)), bottom, bottom)
    } else if (strand == "-") {
        shape_x = c(start+triangle_width, end, end, start+triangle_width, start)
        shape_y = c(top, top, bottom, bottom, mean(c(top,bottom)))
    } else {
        shape_x = c(start, end, end, start)
        shape_y = c(top, top, bottom, bottom)
    }
    return(tibble(x=shape_x, y=shape_y, name=name))
}

# pull an attribute field named $field from a data frame that was read using gffRead
# attribute fields occur in the ninth column of a gff file, and are separated
# using $attrsep
getAttributeField = function(x, field, attrsep = ";") {
    s = strsplit(x, split = attrsep, fixed = TRUE)
    sapply(s, function(atts) {
        a = strsplit(atts, split = "=", fixed = TRUE)
        m = match(field, sapply(a, "[", 1))
        if (!is.na(m)) {
            rv = a[[m]][2]
            return(rv)
        }

        a = strsplit(atts, split = " ", fixed = TRUE)
        m = match(field, sapply(a, "[", 1))
        if (!is.na(m)) {
            rv = a[[m]][2]
            return(rv)
        }

        else {
            rv = as.character("unknown")
        }
        return(as.character("unknown"))
    })
}

# read a gff file and return it as a data frame
gffRead = function(gffFile, nrows = -1) {
    cat("Reading ", gffFile, ".  ", sep="")
    gff = read.table(
        gffFile,
        sep="\t",
        as.is=TRUE,
        quote="",
        header=FALSE,
        comment.char="#",
        nrows = nrows,
        colClasses=c("character", "character", "character", "integer","integer",
             "character", "character", "character", "character")
    )
    colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
    cat("Found", nrow(gff), "rows with classes:",
         paste(sapply(gff, class), collapse=", "), "\n")
    return(gff)
}

# plot a given locus' data
plot_locus = function(signal_df, plotStart, plotEnd, chr_name,
                     feats_df = NULL,
                     lineAlpha=1, lineSize=1.25,
                     yvar="log2fc",
                     color_var="strain",
                     color_vals=NULL,
                     linetype_var=NULL,
                     ylabel="",
                     feat_var = "locus_tag",
                     feat_fill_var = "feature",
                     facet = NULL,
                     ylims="detect", plotFeatures=TRUE, arrowLength=1.5,
                     plotMotifLocs=FALSE, motifs_df=NULL,
                     strand_colors=c("#4667FC", "#FB4B13"),
                     feat_types = c("CDS", "tRNA", "rRNA"),
                     feat_colors = c("#43BA37", "#619CFF", "#F0766D"),
                     log_y=FALSE, name_angle=0, upper=NULL, lower=NULL
) {

    yvar = sym(yvar)
    sym_color_var = sym(color_var)
    if (!is.null(linetype_var)) {
        linetype_var = sym(linetype_var)
    }
    print(color_vals)
 
    plot_sig = signal_df %>%
        dplyr::filter(seqname==chr_name, start<=plotEnd, end>=plotStart) %>%
        mutate(position = (start + end) / 2)

    if (plotMotifLocs) {
        plotMotifs = motifs_df %>% 
            dplyr::filter(seqname==chr_name, end > plotStart, start < plotEnd)
    }

    min_sig = min(plot_sig[,yvar])
    max_sig = max(plot_sig[,yvar])
    sig_range = max_sig - min_sig
    geneTop = min_sig - 0.05*sig_range
    geneBottom = geneTop - 0.075*sig_range
    shiftTop = geneBottom - (0.2 * (geneTop - geneBottom))
    shiftBottom = shiftTop - 0.075*sig_range
    geneText = geneBottom - 0.1 * sig_range
    shiftText = shiftBottom - 0.1 * sig_range
    shiftDist = geneBottom - shiftBottom

    plot = ggplot()

    if (plotMotifLocs) {
        #print(plotMotifs %>% head)
        plot = plot +
            geom_rect(
                data=plotMotifs,
                aes(
                    xmin=start/1e6,
                    xmax=end/1e6,
                    ymin=geneTop,
                    ymax=Inf
                ),
                alpha=0.4
            )
    }

    if (!is.null(lower)) {
        if (is.null(upper)) {
            stop("When specifying a lower conf limit, you must also specify an upper conf limit. exiting now.")
        }
    }

    if (!is.null(upper)) {
        upper = sym(upper)
        if (!is.null(lower)) {
            lower = sym(lower)
        } else {
            stop("When specifying an upper conf limit, you must also specify a lower conf limit. exiting now.")
        }
        plot = plot +
            geom_ribbon(
                data = plot_sig,
                aes(
                    x=position/1e6,
                    ymax=!!upper,
                    ymin=!!lower,
                    fill=!!sym_color_var
                ),
                alpha=0.3
            )

        g = ggplot_build(plot)
        fill_vals = unique(g$data[[1]]['fill'])
        if (is.null(color_vals)) {
            plot = plot +
               scale_fill_manual(name=color_var, values=fill_vals$fill) +
               scale_color_manual(name=color_var, values=fill_vals$fill)
        } else {
            distinct_samples = unlist(unique(plot_sig[,color_var]), use.names=FALSE)
            plot = plot +
               scale_fill_manual(name=color_var, breaks=distinct_samples, values=color_vals) +
               scale_color_manual(name=color_var, breaks=distinct_samples, values=color_vals)
        }
    }
   
    if (plotFeatures) {
        if (is.null(feats_df)) {
            stop("You set plotFeatures = TRUE, but did not provide a dataframe of features to plot!")
        }
        plot_feats = feats_df %>% 
            dplyr::filter(seqname==chr_name, end > plotStart, start < plotEnd) %>%
            mutate(midpoint = as.numeric((start + end)/2))
        feat_var = sym(feat_var)
        feat_fill_var = sym(feat_fill_var)
        feat_labs = feat_types

        shift_feats = c("indep", "tRNA")
        non_shift_feats = c("CDS", "5-UTR", "3-UTR", "rRNA")

        names = plot_feats %>% select(!!feat_fill_var) %>% unlist(use.names=F)
        this_nc = names %in% shift_feats

        any_nc = any(this_nc)
        any_non_nc = any(names %in% non_shift_feats)

        shift_ncRNA = FALSE
        if (any_nc && any_non_nc) {
            shift_ncRNA = TRUE
        }

        plot_feats = plot_feats %>%
            mutate(
                is_nc = this_nc,
                text_tmp_y = ifelse(
                    (shift_ncRNA) & (this_nc),
                    shiftText,
                    geneText
                ),
                text_y = ifelse(shift_ncRNA, text_tmp_y-shiftDist, text_tmp_y)
            )


        plot_geoms = plot_feats %>%
            mutate(poly_name = !!feat_var, poly_fill = !!feat_fill_var) %>%
            nest(data=-c(!!feat_var, !!feat_fill_var)) %>%
            mutate(
                fill = !!feat_fill_var,
                shape_data = purrr::map(
                    data,
                    make_polygons,
                    shift_feats = shift_feats,
                    shift_ncRNA = shift_ncRNA,
                    feat_var = feat_var,
                    plotStart = plotStart,
                    plotEnd = plotEnd,
                    shiftTop = shiftTop,
                    shiftBottom = shiftBottom,
                    geneTop = geneTop,
                    geneBottom = geneBottom
                )
            ) %>% select(fill, shape_data) %>%
            unnest(shape_data)

        plot = plot + 
            new_scale_fill() +
            geom_polygon(
                data=plot_geoms,
                aes(
                    x=x/1e6,
                    y=y,
                    group=!!feat_var,
                    fill=fill
                ),
                color="black"
            ) +
            scale_fill_manual(
                values=feat_colors,
                labels=feat_labs,
                breaks=feat_types
            ) +
            geom_text(
                data=plot_feats,
                aes(
                    x=midpoint/1e6,
                    y=text_y,
                    label=!!feat_var
                ),
                parse=T,
                angle=name_angle
            )

        #plot_geoms = plot_feats %>%
        #    rowwise() %>%
        #    summarize(
        #        create_feature_polygon(
        #            !!feat_var, start, end, strand,
        #            top=geneTop, bottom=geneBottom,
        #            plot_width=plotEnd-plotStart
        #        ),
        #        fill = !!feat_fill_var
        #    )

        #print(plot_feats)
        #stop()
        #print(plot_geoms, n=Inf)
        #stop()
        #print(unique(plot_feats$feature))
        #feat_colors = c("#43BA37", "#619CFF", "#F0766D")
        #feat_types = c("CDS", "tRNA", "rRNA")

        #plot = plot + 
        #    new_scale_fill() +
        #    geom_polygon(
        #        data=plot_geoms,
        #        aes(
        #            x=x/1e6,
        #            y=y,
        #            group=!!feat_var,
        #            fill=fill
        #        ),
        #        color="black"
        #    ) +
        #    #scale_fill_manual(values=strand_colors) +
        #    geom_text(
        #        data=plot_feats,
        #        aes(
        #            x=midpoint/1e6,
        #            y=geneText,
        #            label=!!feat_var
        #        ),
        #        parse=T,
        #        angle=name_angle
        #    )
        #if (!is.null(feat_colors)) {
        #    plot = plot + scale_fill_manual(values=feat_colors, labels=feat_types, breaks=feat_types)
        #}
    }

    #print(plot_sig %>% head)
    plot = plot + 
        geom_line(
            data=plot_sig,
            aes(
                x=position/1e6,
                y=!!yvar,
                color=!!sym_color_var,
                linetype=!!linetype_var
            ),
            size=lineSize,
            alpha=lineAlpha
        ) +
        #scale_color_manual(values=sample_colors) +
        theme_classic() +
        theme(
            text = element_text(size=10),
            axis.line = element_line(size=1),
            axis.text = element_text(size=9, color="black"),
            axis.ticks = element_line(color="black")
        ) +
        labs(y=ylabel, x="Genome position (Mb)") +
        coord_cartesian(xlim = c(plotStart/1e6, plotEnd/1e6)) +
        guides(colour = guide_legend(override.aes = list(alpha = 1)))

        if (is.null(color_vals)) {
            plot = plot +
               scale_color_manual(name=color_var, values=fill_vals$fill)
        } else {
            distinct_samples = unlist(unique(plot_sig[,color_var]), use.names=FALSE)
            plot = plot +
               scale_color_manual(name=color_var, breaks=distinct_samples, values=color_vals)
        }

    if (!(ylims == "detect")) {
        plot = plot + coord_cartesian(
            xlim = c(plotStart/1e6, plotEnd/1e6),
            ylim=ylims
        )
    }

    if (!is.null(facet)) {
        facets = str_split(str_remove_all(facet, " "), "~", simplify=TRUE)
        row_var = facets[1]
        col_var = facets[2]
        if (row_var == "") {
            row_var = NULL
        } else {
            print(row_var)
            row_var = sym(row_var)
            row_var = vars(!!row_var)
        }
        if (col_var == "") {
            col_var = NULL
        } else {
            print(col_var)
            col_var = sym(col_var)
            col_var = vars(!!col_var)
        }
        plot = plot + facet_grid(rows=row_var, cols=col_var)
        #facet = as.formula(facet)
        #plot = plot + facet_grid(facet)
    }
    if (log_y) {
        plot = plot + scale_y_log10()
    }

    return(plot)
}


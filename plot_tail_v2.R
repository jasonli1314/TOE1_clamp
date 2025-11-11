#Some old scripts or packages set the RNG kind to "Rounding" for compatibility
RNGkind(sample.kind="Rejection")
library("ggplot2")
library("argparse")
#library("stringr")
#library("ggseqlogo")

get_flag_input <- function() {
	parser <- ArgumentParser()
	parser$add_argument("-l", "--legend", action="store_true",
	                    default=F, help="Show legend [default off]")

	parser$add_argument("-s", "--graph_size", type="double",
	                    default=c(3, 2), nargs=2, metavar=c('w','h'),
	                    help="Height and width of graph [default: c(3, 2)]")

	parser$add_argument("-x", "--x_range", type="integer",
	                    default=c(-25, 10), nargs=2, metavar=c('x-min','x-max'),
	                    help="Plot x range [default: -5 10]")
	
	parser$add_argument("-y", "--y_range", type="integer",
	                    default=c(0, 100), nargs=2, metavar=c('y-min','y-max'),
	                    help="Plot y range [default: 0 100]")
		
	parser$add_argument("-a", "--analysis_range", type="integer",
	                    default=c(-50, 50), nargs=2, metavar=c('min','max'),
	                    help="Analysis range (should be >= plot range) [default: -50 50]")
	
	parser$add_argument("-i", "--in_folder", type="character",
	                    default="", help="input file folder")
	
	parser$add_argument("-o", "--out_folder", type="character",
	                    default="", help="Output file folder")
	
	parser$add_argument("-e", "--experiment_manifest",
	                    type="character", default="",
	                    help="Path to experiment file manifest")
	
	parser$add_argument("-g", "--gene_manifest",
	                    type="character", default="",
	                    help="Path to gene information manifest")
	
	args <- parser$parse_args()
	return(args)
}


get_exp_colors <- function(experimentinfo_df) {
  exps <- unique(experimentinfo_df$expname)
  color_vec <- experimentinfo_df$color[match(exps, experimentinfo_df$expname)]
  names(color_vec) <- exps
	return(color_vec)
}


upload_data <- function(experimentinfo_df, RNA) {
  experimentinfo_df <- subset(experimentinfo_df, RNA_name==RNA)
	data_list <- list()
	#Needed to skip the header; missing value
	#Specifying quote is important for not stumbling on 3'end comments
	for (i in 1:nrow(experimentinfo_df)) {
	  input_file <- file.path(IN_FOLDER, experimentinfo_df$tailfile[i])
		df <- read.table(input_file, header=F, sep=",", quote="\"", skip=1, fill=T, stringsAsFactors=F) 
		names(df) <- c("Sequence", "Unique_Reads", "Gene", "3_end", "Tail_length", "Tail_seq", "Notes")
		data_list <- c(data_list, list(df))}
	experimentinfo_df$data <- data_list
	return(experimentinfo_df)
}


extract_gene_specific_taildata <- function(taildata_df, plotgeneids, end_correction) {
  
	#Extracts lines from taildata files matching current gene
	filtered_taildata_df <- data.frame(matrix(NA,0,7), stringsAsFactors=F)
	names(filtered_taildata_df) <- c("Sequence", "Unique_Reads", "Gene", "3_end", "Tail_length", "Tail_seq", "Notes")

	for (n in 1:nrow(taildata_df)) {
		read_info <- taildata_df[n,]
		temp_line <- as.character(read_info$Gene)
		#Placing all ENS numbers in the genes into a vector
		geneids <- regmatches(temp_line, gregexpr("ENS.[0-9]{11}", temp_line))[[1]]	
		#Adding the row to the plotgene_df if all geneid numbers in the specific data line 
		# are in the geneids to be plotted (plotgeneids) and the 3' end is in the analysis range
		if (length(geneids) == 0) { next }
		if (all(geneids %in% plotgeneids)) {
		  #This corrects for any misannotated 3' end
			read_info$"3_end" <- as.numeric(read_info$"3_end") - end_correction
			if (read_info$"3_end" >= ANALYSIS_RANGE[1] & read_info$"3_end" <= ANALYSIS_RANGE[2]) {
				filtered_taildata_df <- rbind(filtered_taildata_df, read_info, stringsAsFactors=F)
			}
		}
	}

	return(filtered_taildata_df)
}


calculate_endpos <- function(temp_df) {

	end_df <- data.frame(short=c(rep(0, ANALYSIS_WINDOW)), full=c(rep(0, ANALYSIS_WINDOW)))

	for (n in 1:nrow(temp_df)) {
		line_reads <- temp_df$Unique_Reads[n]
		#Calculating end positions
		short_end_pos <- temp_df$"3_end"[n]
		full_end_pos <- temp_df$"3_end"[n] + temp_df$Tail_length[n]
		full_end_pos <- min(full_end_pos, ANALYSIS_RANGE[2])
		
		#Calculating row positions using adjustment factor
		short_row_pos <- short_end_pos + POS_ADJ
		full_row_pos <- full_end_pos + POS_ADJ
		end_df$short[short_row_pos] <- end_df$short[short_row_pos] + line_reads
		end_df$full[full_row_pos] <- end_df$full[full_row_pos] + line_reads	
	}
	return(end_df)
}


get_tail_comp <- function(temp_df, LOGO_MAX_POS=8) {
  #Setting up a temp matrix to hold the PT-tail data for each experiment
  tail_m <- matrix(0, 4, 8)
  rownames(tail_m) <- c("A", "T", "C", "G")
  
  for (row in 1:nrow(temp_df)) {
    reads <- temp_df$Unique_Reads[row]
    tail_length <- temp_df$Tail_length[row]
    if (tail_length > LOGO_MAX_POS) {tail_length <- LOGO_MAX_POS}
    
    if (tail_length > 0) {
      for (n in 1:tail_length) {
        nuc <- substr(temp_df$Tail_seq[row], n, n)
        if (nuc %in% rownames(tail_m)) { #needed to evade N
          tail_m[nuc, n] <- tail_m[nuc, n] + reads
        }
      }
    }
  }
  rownames(tail_m) <- c("A", "U", "C", "G")
  return(tail_m)
}

CDF_plot <- function(plot_df, linecolors, filename) {

	#xminorticks <- seq(PLOT_RANGE[1], PLOT_RANGE[2], by=5)
	yminorticks <- seq(PCT_RANGE[1], PCT_RANGE[2], by=20)
	
	p <- ggplot(plot_df) + coord_cartesian() + theme_classic() + ggtitle(sub(".pdf", "", basename(filename))) + 

  annotate("rect", xmin=PLOT_RANGE[1], xmax=0, ymin=PCT_RANGE[1], ymax=PCT_RANGE[2], fill="grey", alpha=0.2) + 

	scale_x_continuous(name="3' end position", expand=c(0, 0), limits=c(PLOT_RANGE[1], PLOT_RANGE[2])) +

	scale_y_continuous(name="Cumulative Percent", breaks=yminorticks, 
		expand=c(0, 0), limits=c(PCT_RANGE[1], PCT_RANGE[2]), position="left") +

	geom_step(aes(x=Position, y=Percent_Full, color=Experiment, group=Experiment), linewidth=0.7, 
	  alpha=0.9, show.legend=LEGEND) + 
	
	geom_step(aes(x=Position, y=Percent_Short, color=Experiment, group=Experiment), linewidth=0.7, linetype=3, 
		alpha=0.9, show.legend=F) +
	 
	scale_color_manual(values=linecolors) +

	geom_rect(aes(xmin=Position, xmax=Position+0.99, ymin=Percent_Full, ymax=Percent_Short, color=Experiment, 
	              group=Experiment, fill=Experiment), linetype=0, alpha=0.2, show.legend=F) +

	scale_fill_manual(values=linecolors)
	
	ggsave(filename, plot=p, width=WIDTH, height=HEIGHT)

}

get_plot_data <- function(experimentinfo_df, plotgeneids, end_correction){
  CDF_plot_df <- data.frame(matrix(NA,0,4), stringsAsFactors=F)
  names(CDF_plot_df) <- c("Experiment", "Position", "Percent_short", "Percent_full")
  PMF_plot_df <- CDF_plot_df
  logo_m <- matrix(NA, 0, LOGO_MAX_POS)
  raw_list <- list()
  res_list <- list()

  #looping over each experiment/condition (including a few replicates) per plot
  for (exp_name in unique(experimentinfo_df$expname)) {
    
    curr_exp_df <- subset(experimentinfo_df, expname==exp_name)
    compiled_shortend_df <- data.frame(matrix(NA, ANALYSIS_WINDOW, 0), stringsAsFactors=F)
    compiled_fullend_df <- data.frame(matrix(NA, ANALYSIS_WINDOW, 0), stringsAsFactors=F)
    temp_logo_m <- 0
    
    for (f in 1:nrow(curr_exp_df)) {
      temp_df <- extract_gene_specific_taildata(curr_exp_df$data[[f]], plotgeneids, end_correction)
      # Update raw data per replicate
      raw_list <- c(raw_list, list(temp_df))
      
      if (nrow(temp_df) >	0) {
        total_reads <- sum(temp_df$Unique_Reads)
        # End position CDF data
        end_df <- calculate_endpos(temp_df)
        compiled_shortend_df <- cbind(compiled_shortend_df, end_df$short/total_reads*100)
        compiled_fullend_df <- cbind(compiled_fullend_df, end_df$full/total_reads*100)
        # Tail Logo data
        temp_logo_m <- temp_logo_m + get_tail_comp(temp_df)/total_reads*100
      } else { next }
    }
    
    #Update PMF/CDF_plot_df for each exp (average of ~3 biological replicates)
    cum_shortend_df <- data.frame(percent=cumsum(rowMeans(compiled_shortend_df)))
    cum_fullend_df <- data.frame(percent=cumsum(rowMeans(compiled_fullend_df)))
    temp_CDF_df <- data.frame("Experiment"=c(rep(exp_name, (PLOT_RANGE[2]-PLOT_RANGE[1]+1))), 
                               "Position"=c(PLOT_RANGE[1]:PLOT_RANGE[2]), 
                               "Percent_Short"=cum_shortend_df$percent[(PLOT_RANGE[1]+POS_ADJ):(PLOT_RANGE[2]+POS_ADJ)], 
                               "Percent_Full"=cum_fullend_df$percent[(PLOT_RANGE[1]+POS_ADJ):(PLOT_RANGE[2]+POS_ADJ)])
    
    
    mean_shortend_df <- data.frame(percent=rowMeans(compiled_shortend_df))
    mean_fullend_df <- data.frame(percent=rowMeans(compiled_fullend_df))
    temp_PMF_df <- data.frame("Experiment"=c(rep(exp_name, (PLOT_RANGE[2]-PLOT_RANGE[1]+1))), 
                              "Position"=c(PLOT_RANGE[1]:PLOT_RANGE[2]), 
                              "Percent_Short"=mean_shortend_df$percent[(PLOT_RANGE[1]+POS_ADJ):(PLOT_RANGE[2]+POS_ADJ)], 
                              "Percent_Full"=mean_fullend_df$percent[(PLOT_RANGE[1]+POS_ADJ):(PLOT_RANGE[2]+POS_ADJ)])
    
    CDF_plot_df <- rbind(CDF_plot_df, temp_CDF_df)
    PMF_plot_df <- rbind(PMF_plot_df, temp_PMF_df)
    
    # Update logo_m for each exp
    logo_m <- rbind(logo_m, temp_logo_m/nrow(curr_exp_df))
  }
  
  res_list$CDF_plot_df <- CDF_plot_df
  res_list$PMF_plot_df <- PMF_plot_df
  
  res_list$logo_m <- logo_m
  res_list$raw_list <- raw_list
  return(res_list)
}

###########################################
#Unpacking analysis range input:
args <- get_flag_input()
LOGO_MAX_POS <- 8
ANALYSIS_RANGE <- args$analysis_range
#Adjustment value to convert 3' end position to row number in analysis data frames
POS_ADJ <- 1 - (ANALYSIS_RANGE[1])
ANALYSIS_WINDOW <- ANALYSIS_RANGE[2] - ANALYSIS_RANGE[1] + 1
LEGEND <- args$legend
WIDTH <- args$graph_size[1]
HEIGHT <- args$graph_size[2]
PLOT_RANGE <- args$x_range
PCT_RANGE <- args$y_range
IN_FOLDER <- args$in_folder
OUT_FOLDER <- args$out_folder
GENE_PATH <- args$gene_manifest
EXP_PATH <- args$experiment_manifest

experimentinfo_raw <- read.csv(EXP_PATH)
RNA_vec <- unique(experimentinfo_raw$RNA_name)

# Loop over different RNA types
for (RNA in RNA_vec){
  
  #Reading datafiles into the experimentinfo data frame
  experimentinfo_df <- upload_data(experimentinfo_raw, RNA)
  linecolors <- get_exp_colors(experimentinfo_df)
  geneinfo_df <- read.csv(GENE_PATH, row.names=1)
  #extracting geneid(s) to be plotted into a vector (plotgeneids)
  #(there could be multiple since some RNAs are from repeat genes)
  plotgeneids <- strsplit(geneinfo_df[RNA, "geneid"], ", ")[[1]]
  #extracting the 3' end position [in case it is misannotated]
  end_correction <- geneinfo_df[RNA, "endpos"]
  
  res_list <- get_plot_data(experimentinfo_df, plotgeneids, end_correction)
  CDF_plot_df <- res_list$CDF_plot_df
  PMF_plot_df <- res_list$PMF_plot_df
  logo_m <- res_list$logo_m
  raw_list <- res_list$raw_list
  
  if (!dir.exists(file.path(OUT_FOLDER, 'CDF'))) {dir.create(file.path(OUT_FOLDER, 'CDF'), recursive=T)}
  if (!dir.exists(file.path(OUT_FOLDER, 'PMF'))) {dir.create(file.path(OUT_FOLDER, 'PMF'), recursive=T)}
  if (!dir.exists(file.path(OUT_FOLDER, 'logo'))) {dir.create(file.path(OUT_FOLDER, 'logo'), recursive=T)}
  if (!dir.exists(file.path(OUT_FOLDER, 'raw'))) {dir.create(file.path(OUT_FOLDER, 'raw'), recursive=T)}
  
  last_exp <- tail(experimentinfo_df$expname, 1)
  CDF_filename <- file.path(OUT_FOLDER, 'CDF', paste(last_exp, RNA, "CDF.pdf", sep="_"))
  PMF_filename <- file.path(OUT_FOLDER, 'PMF', paste(last_exp, RNA, "PMF.Rdat", sep="_"))
  logo_filename <- file.path(OUT_FOLDER, 'logo', paste(last_exp, RNA, "logo.csv", sep="_"))
  raw_filename <- file.path(OUT_FOLDER, 'raw', paste(last_exp, RNA, "raw.Rdat", sep="_"))
  
  CDF_plot(CDF_plot_df, linecolors, CDF_filename)
  write.csv(logo_m, logo_filename)
  save(PMF_plot_df, file=PMF_filename)
  save(raw_list, file=raw_filename)
  
}





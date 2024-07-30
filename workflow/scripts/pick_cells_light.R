library(argparse)
library(ggplot2)
library(patchwork)
library(caret)
library(pheatmap)
library(reshape2)
library(ComplexUpset)
library(IRanges)

# Parse command line arguments
parser <- ArgumentParser()
parser$add_argument("--metrics_all", help="Input file", required=TRUE)
parser$add_argument("--metrics_peak", help="Input file", required=TRUE)
parser$add_argument("--out", help="Output metadata table with picked cells", required=TRUE)
parser$add_argument("--cells_list", help="Output file with cell barcodes of picked cells", required=TRUE)
parser$add_argument("--kmeans", help="Number of clusters for kmeans", default=2)
parser$add_argument("--knee_stdev", help="Number of standard deviations above the mean for knee detection", default=3, type='numeric')
parser$add_argument("--skip_first", help="Number of cells with the highest numbers of reads to skip during sd thresholding [typically doublets/clumps]\n Can be eihter percentile [0-0.99] or absolute number [e.g. 500]", default=0.1, type='numeric')
parser$add_argument("--minimum_cells", help="Minimum number of cells to pick", default=1000, type='numeric')
args <- parser$parse_args()

out_prefix <- dirname(args$out)


# Functions 
find_intervals_in_vector <- function(x){

#    Takes vector as input e.g. '1''2''3''4''5''6''7''8''3379''3380''3381''3382''3383''3384'

#    Returns interval starts and ends:
#    interval_starts	interval_ends
#    1	                8
#    3379	            3384
    
    # If there is a knee below 100 cells, remove it 
    x <- x[x > 100]
    if(length(x) == 0){
        return(IRanges())
    }
    
    # Merge intervals if difference between them is less than 50 cells
    interval_starts = c(x[x==min(x)], x[which(diff(x) > 1) + 1])
    interval_ends   = c(x[which(diff(x) > 1)],max(x))
    x.iranges <- IRanges::IRanges(start=interval_starts,end=interval_ends)
    x.iranges <- IRanges::reduce(x.iranges, min.gapwidth=50)
    return(as.data.frame(x.iranges))
}


find_clusters_kmeans <- function(metrics.df, nclusters = 2){
    # Find the clusters using kmeans
    metrics.df$cluster <- as.factor(kmeans(as.matrix(metrics.df[,c('FRiP','log_all_reads')]),centers = nclusters)$cluster)
    # Rename the clusters
    log_all_summary           <- aggregate(metrics.df$log_all_reads,by=list(metrics.df$cluster),FUN=median)
    log_all_max               <- which(log_all_summary$x == max(log_all_summary$x))
    metrics.df$is_cell_kmeans <- FALSE
    metrics.df$is_cell_kmeans[metrics.df$cluster == log_all_max] = TRUE
    return(metrics.df)
}

find_clusters_knee_plot <- function(metrics.df, knee_stdev = 3){
    metrics.df               <- metrics.df[order(metrics.df$all_reads,decreasing = TRUE),]

    if(args$skip_first < 1){
        i.skip <- which(metrics.df$all_reads > quantile(metrics.df$all_reads,1-args$skip_first)) 
    }
    else{
        i.skip <- 1:args$skip_first
    }
    
    # Calculate diff between ranked cells and smoothen the diff to remove noise/outliers
    reads.diff.smooth.full        <- runmed(-diff(metrics.df$log_all_reads),k=101)
    reads.diff.smooth.truncated   <- runmed(-diff(metrics.df$log_all_reads[-i.skip]),k=101) # Truncated 
 
    # Find peaks in the smoothened diff plot
    peaks <- which(reads.diff.smooth.truncated/sd(reads.diff.smooth.truncated) > args$knee_stdev)
    peak_intervals = find_intervals_in_vector(peaks)

    
    if( is.null(nrow(peak_intervals)) ){
        warning('INFO: No knees found using the knee plot method. Consider increasing the stdev threshold for knee detection')
        warning('INFO: Manualy setting knee to 5000 cells')
        knee      = 5000
        all_knees = knee
    }
    else if(nrow(peak_intervals) > 1){
        warning('INFO: More than one knee found using the knee plot method. Using the first knee as the cutoff\nConsider increasing the stdev threshold for knee detection')
        all_knees = peak_intervals$start
    }
    else if (nrow(peak_intervals) == 1){
        all_knees = peak_intervals$start
    }
    all_knees = all_knees + max(c(0,i.skip))
    knee      = all_knees[1]    
    
    # Take the first knee as the cutoff
    metrics.df$is_cell_knee = FALSE
    metrics.df$is_cell_knee[1:knee] = TRUE
    
    # Plot some QC
    # Recalculate the diff and smoothened diff for plotting
    
    plot.df      <- data.frame(sd = reads.diff.smooth.full / sd(reads.diff.smooth.truncated), rank=1:length(reads.diff.smooth.full))
    p1 <- ggplot(data=plot.df) + 
      geom_point(aes(x=rank,y=sd)) + 
      geom_hline(yintercept = args$knee_stdev,color='red') + 
      geom_vline(xintercept = all_knees,color='blue') +
      geom_vline(xintercept = knee,color='darkgreen',lwd=2,linetype='dashed') +
      ggtitle(paste0('Smoothened diff of ranked cells / sd\n', 
                      '  Red horizontal line = stdev threshold for knee calling (adjust if needed using --knee_stdev)\n',
                      '  Blue vertical lines = detected knees\n',
                      '  Green dashed vertical line = used knee cutoff\n',
                      '  ncells=', sum(metrics.df$is_cell_knee),'\n')) + 
      xlab('Cell barcode rank') + ylab('Smoothened diff in read number / sd [zscore]')
    ggsave(paste0(out_prefix,'/knee_plot_QC.pdf'),plot=p1,width=10,height=8)


    return(metrics.df)
}

# Load the data 
metrics.df <- merge(read.table(file=args$metrics_all), read.table(file=args$metrics_peak),by='V1')
colnames(metrics.df) <- c('barcode','all_reads','peak_reads')

# Filter away droplets with < 50 reads
metrics.df <- metrics.df[metrics.df$all_reads > 50,]

# Add cell rank
metrics.df               <- metrics.df[order(metrics.df$all_reads,decreasing = TRUE),]
metrics.df$rank          <- 1:nrow(metrics.df)

# Add FRiP
metrics.df$FRiP          <- metrics.df$peak_reads/metrics.df$all_reads

# Log Transform all reads for convenience 
metrics.df$log_all_reads <- log10(metrics.df$all_reads)

################################
# Find using kmeans clustering #
################################

metrics.df <- find_clusters_kmeans(metrics.df, nclusters = args$kmeans)

########################
# Find using knee plot #
########################

metrics.df <- find_clusters_knee_plot(metrics.df, knee_stdev = args$knee_stdev)


# Save the QC plots

options(repr.plot.width=20,repr.plot.height=8)
p1 <- ggplot(data=metrics.df,aes(x=all_reads,y=peak_reads/all_reads)) + 
  geom_point(aes(col=is_cell_kmeans),size=0.5) + 
  scale_x_log10() + ggtitle(paste0('Cells picking using kmeans clustering\n  Ncells=',sum(metrics.df$is_cell_kmeans))) 

p2 <- ggplot(data=metrics.df,aes(x=all_reads,y=peak_reads/all_reads)) + 
  geom_point(aes(col=is_cell_knee),size=0.5) + 
  scale_x_log10() + 
  ggtitle(paste0('Cell picking using knee plot\n  Ncells=',sum(metrics.df$is_cell_knee)))
p3 <- ggplot(data=metrics.df,aes(x=as.numeric(rank),y=all_reads)) + 
  geom_point(aes(color=is_cell_knee)) + 
  scale_y_log10() + 
  ylab('Number of reads per cell') + xlab('Cell barcode rank')
  

# Confusion matrix and upset plot between the 2 methods
cm      <- confusionMatrix(as.factor(metrics.df$is_cell_knee),as.factor(metrics.df$is_cell_kmeans)) 
cm.long <- reshape2::melt(cm$table)
colnames(cm.long) <- c('is_cell_knee','is_cell_kmeans','value')
print(cm.long)

p4 <- upset(metrics.df[,c('is_cell_knee','is_cell_kmeans')],intersect=c('is_cell_knee','is_cell_kmeans'))
p5 <- ggplot(data=cm.long,aes(x=is_cell_knee,y=is_cell_kmeans)) + geom_tile(aes(fill=value)) + geom_text(aes(label = value)) + scale_fill_gradient(low='white',high='red')


pdf(paste0(out_prefix,'/Cell_picking_QC_plot.pdf'), width=10, height=8)
print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
dev.off()
rownames(metrics.df) <- metrics.df$barcode

# Export the metadata table
write.csv(file=args$out,x=metrics.df)

# Export the list of picked cells 
write.table(file = args$cells_list, 
            x = metrics.df$barcode[metrics.df$is_cell_knee & metrics.df$is_cell_kmeans], 
            row.names = FALSE, col.names = FALSE, 
            quote = FALSE)

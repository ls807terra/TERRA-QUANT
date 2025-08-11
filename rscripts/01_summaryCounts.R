#!/staging/biology/ls807terra/0_Programs/anaconda3/envs/RNAseq_quantTERRA/bin/Rscript

## Options
pacman::p_load("optparse")

option_list = list(
  make_option(c("-c", "--counts"), type="character", default=NULL,
              help="Enter a directory that contains count files.", metavar="COUNTS"),
  make_option(c("-o", "--output"), type="character", default="./raw_counts_table",
              help="Specify a output count table file. [default= %default]", metavar="OUTPUT"),
  make_option(c("-r", "--telo_repeat"), type="character", default="./TERRA_repeat_table",
              help="Specify a output count table file. [default= %default]", metavar="OUTPUT"),
  make_option(c("-s", "--subtelo"), type="character", default="./TERRA_subtelo_table",
              help="Specify a output count table file. [default= %default]", metavar="OUTPUT"),            
  make_option(c("-f", "--format"), type="character", default="csv",
              help="Specify output file format. Possible: csv xlsx [default= %default]", metavar="FORMAT"),
  make_option(c("-k", "--keyword"), type="character", default=".count",
              help="If the count files did not name as *.count*, specify a common name of those count files. 
		    Example: .txt [default= %default]", metavar="KEYWORD")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Option validator
if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("Please Enter at least -c option: -c count_path/", call.=FALSE)
}
if (opt$format != "csv" & opt$format != "xlsx"){
  stop("Output format only support csv or xlsx", call.=FALSE)
}

## Loading packages
pacman::p_load(rio)

## Loading files
countDir <- opt$counts
files <- paste0(countDir,"/",grep(opt$keyword,list.files(countDir),value = T))

geneNum <- 60945

counts.df <- import(files[1],format = "csv")[1:geneNum,]
col.index <- 2
for (i in files) {
  tmp.df <- import(i,format = "csv")[1:geneNum,]
  counts.df[,col.index] <- tmp.df[,2]
  col.index <- col.index + 1
  print(paste("Process",files,"...",col.index))
}
colnames(counts.df) <- c("Gene",grep(".count",list.files(countDir),value = T))

genes_to_exclude <- c(
  "TERRA-chr10p-repeat", "TERRA-chr10p-subtelo", "TERRA-chr10q-repeat", "TERRA-chr10q-subtelo",
  "TERRA-chr11p-repeat", "TERRA-chr11p-subtelo", "TERRA-chr11q-repeat", "TERRA-chr11q-subtelo",
  "TERRA-chr12p-repeat", "TERRA-chr12p-subtelo", "TERRA-chr13q-repeat", "TERRA-chr13q-subtelo",
  "TERRA-chr14q-repeat", "TERRA-chr14q-subtelo", "TERRA-chr15p-repeat", "TERRA-chr15p-subtelo",
  "TERRA-chr15q-repeat", "TERRA-chr15q-subtelo", "TERRA-chr16p-repeat", "TERRA-chr16p-subtelo",
  "TERRA-chr16q-repeat", "TERRA-chr16q-subtelo", "TERRA-chr17p-repeat", "TERRA-chr17p-subtelo",
  "TERRA-chr17q-repeat", "TERRA-chr17q-subtelo", "TERRA-chr18p-repeat", "TERRA-chr18p-subtelo",
  "TERRA-chr18q-repeat", "TERRA-chr18q-subtelo", "TERRA-chr19p-repeat", "TERRA-chr19p-subtelo",
  "TERRA-chr19q-repeat", "TERRA-chr19q-subtelo", "TERRA-chr1p-repeat", "TERRA-chr1p-subtelo",
  "TERRA-chr1q-repeat", "TERRA-chr1q-subtelo", "TERRA-chr20p-2-subtelo", "TERRA-chr20p-repeat",
  "TERRA-chr20p-subtelo", "TERRA-chr20q-repeat", "TERRA-chr20q-subtelo", "TERRA-chr21q-repeat",
  "TERRA-chr21q-subtelo", "TERRA-chr22q-subtelo", "TERRA-chr2p-repeat", "TERRA-chr2p-subtelo",
  "TERRA-chr2q-repeat", "TERRA-chr2q-subtelo", "TERRA-chr3p-repeat", "TERRA-chr3p-subtelo",
  "TERRA-chr3q-repeat", "TERRA-chr3q-subtelo", "TERRA-chr4p-repeat", "TERRA-chr4p-subtelo",
  "TERRA-chr4q-repeat", "TERRA-chr4q-subtelo", "TERRA-chr5p-repeat", "TERRA-chr5p-subtelo",
  "TERRA-chr5q-repeat", "TERRA-chr5q-subtelo", "TERRA-chr6p-repeat", "TERRA-chr6p-subtelo",
  "TERRA-chr6q-repeat", "TERRA-chr6q-subtelo", "TERRA-chr7p-repeat", "TERRA-chr7p-subtelo",
  "TERRA-chr7q-repeat", "TERRA-chr7q-subtelo", "TERRA-chr8p-repeat", "TERRA-chr8p-subtelo",
  "TERRA-chr9p-repeat", "TERRA-chr9p-subtelo", "TERRA-chr9q-repeat", "TERRA-chr9q-subtelo",
  "TERRA-chrXq-repeat", "TERRA-chrXq-subtelo", "TERRA_chr2_5_CTAACCn_r2_repeat",
  "TERRA_chr2_CTAACCn_3_r2_repeat", "chr11_ITS", "chr13_ITS", "chr14_ITS_1", "chr14_ITS_2",
  "chr14_ITS_3", "chr15_ITS", "chr18p-2", "chr18p_1", "chr18p_2", "chr19_ITS", "chr1_ITS",
  "chr20_ITS_1", "chr20_ITS_2", "chr22_ITS_1", "chr22_ITS_2", "chr2_ITS_1", "chr2_ITS_2",
  "chr3_ITS", "chr5p-2", "chr7_ITS", "chr8_ITS"
)
telomeric_repeat_counts <- counts.df %>% filter(grepl("repeat", Gene))
subtelomeric_TERRA_counts <- counts.df %>% filter((Gene %in% genes_to_exclude) & !grepl("repeat", Gene))
filtered_counts <- counts.df %>% filter(!(Gene %in% genes_to_exclude))
remove_index <- grep("__",filtered_counts$Gene)
select_i <- seq(0,ncol(filtered_counts),by = 2)
filtered_counts <- filtered_counts[-remove_index,c(1,select_i)]

export(filtered_counts,opt$output,format = opt$format)
export(counts.df,opt$telo_repeat,format = opt$format)
export(subtelomeric_TERRA_counts,opt$subtelo,format = opt$format)


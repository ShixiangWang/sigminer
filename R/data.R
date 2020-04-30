## how to document datasets: you need to specify @docType and @name; do not
## forget NULL in the end

#' Location of Centromeres at Genome Build hg19
#' @docType data
#' @name centromeres.hg19
#' @format A data.frame
#' @source Generate from UCSC gold path
#' @examples
#' data(centromeres.hg19)
NULL

#' Location of Centromeres at Genome Build hg38
#' @docType data
#' @name centromeres.hg38
#' @format A data.frame
#' @source Generate from Genome Reference Consortium
#' @examples
#' data(centromeres.hg38)
NULL

#' Chromosome Size of Genome Build hg19
#' @docType data
#' @name chromsize.hg19
#' @format A data.frame
#' @source Generate from UCSC gold path
#' @examples
#' data(chromsize.hg19)
NULL

#' Chromosome Size of Genome Build hg38
#' @docType data
#' @name chromsize.hg38
#' @format A data.frame
#' @source Generate from UCSC gold path
#' @examples
#' data(chromsize.hg38)
NULL

#' Location of Chromosome Cytobands at Genome Build hg19
#' @docType data
#' @name cytobands.hg19
#' @format A data.frame
#' @source from UCSC
#' @examples
#' data(cytobands.hg19)
NULL

#' Location of Chromosome Cytobands at Genome Build hg38
#' @docType data
#' @name cytobands.hg38
#' @format A data.frame
#' @source from UCSC
#' @examples
#' data(cytobands.hg38)
NULL

#' Merged Transcript Location at Genome Build hg19
#' @docType data
#' @name transcript.hg19
#' @format A `data.table`
#' @source from GENCODE release v33.
#' @examples
#' data(transcript.hg19)
NULL

#' Merged Transcript Location at Genome Build hg38
#' @docType data
#' @name transcript.hg38
#' @format A `data.table`
#' @source from GENCODE release v33.
#' @examples
#' data(transcript.hg38)
NULL

#' Classfication Table of Copy Number Features Devised by Wang et al.
#' @docType data
#' @name CN.features
#' @format A `data.table` with "sigminer.features" class name
#' @source Generate from code under data_raw/
#' @examples
#' data(CN.features)
NULL

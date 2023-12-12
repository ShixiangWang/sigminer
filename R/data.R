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

#' Location of Centromeres at Genome Build T2T
#' @docType data
#' @name centromeres.T2T
#' @format A data.frame
#' @source from T2T study
#' @examples
#' data(centromeres.T2T)
NULL

#' Location of Centromeres at Genome Build mm10
#' @docType data
#' @name centromeres.mm10
#' @format A data.frame
#' @source Generate from <https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/gap.txt.gz>
#' @examples
#' data(centromeres.mm10)
NULL

#' Location of Centromeres at Genome Build mm9
#' @docType data
#' @name centromeres.mm9
#' @format A data.frame
#' @source Generate from <https://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/>
#' with code:
#' ```sh
#' for i in $(seq 1 19) X Y;
#' do
#' wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/chr${i}_gap.txt.gz
#' done
#' ```
#' @examples
#' data(centromeres.mm9)
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

#' Chromosome Size of Genome Build T2T
#' @docType data
#' @name chromsize.T2T
#' @format A data.frame
#' @source from T2T study
#' @examples
#' data(chromsize.T2T)
NULL

#' Chromosome Size of Genome Build mm10
#' @docType data
#' @name chromsize.mm10
#' @format A data.frame
#' @source Generate from UCSC gold path <http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes>
#' @examples
#' data(chromsize.mm10)
NULL
#' Chromosome Size of Genome Build mm9
#' @docType data
#' @name chromsize.mm9
#' @format A data.frame
#' @source Generate from UCSC gold path <http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes>
#' @examples
#' data(chromsize.mm9)
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

#' Location of Chromosome Cytobands at Genome Build T2T
#' @docType data
#' @name cytobands.T2T
#' @format A data.frame
#' @source from T2T study
#' @examples
#' data(cytobands.T2T)
NULL

#' Location of Chromosome Cytobands at Genome Build mm10
#' @docType data
#' @name cytobands.mm10
#' @format A data.frame
#' @source from UCSC <http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/cytoBand.txt.gz>
#' @examples
#' data(cytobands.mm10)
NULL
#' Location of Chromosome Cytobands at Genome Build mm9
#' @docType data
#' @name cytobands.mm9
#' @format A data.frame
#' @source from UCSC <http://hgdownload.cse.ucsc.edu/goldenpath/mm9/database/cytoBand.txt.gz>
#' @examples
#' data(cytobands.mm9)
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

#' Merged Transcript Location at Genome Build T2T
#' @docType data
#' @name transcript.T2T
#' @format A `data.table`
#' @source from T2T study.
#' @examples
#' data(transcript.T2T)
NULL

#' Merged Transcript Location at Genome Build mm10
#' @docType data
#' @name transcript.mm10
#' @format A `data.table`
#' @source from GENCODE release M25.
#' @examples
#' data(transcript.mm10)
NULL
#' Merged Transcript Location at Genome Build mm9
#' @docType data
#' @name transcript.mm9
#' @format A `data.table`
#' @source from UCSC <http://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/transcriptome.txt.gz>
#' @examples
#' data(transcript.mm9)
NULL
#' Classification Table of Copy Number Features Devised by Wang et al. for Method 'W'
#' @docType data
#' @name CN.features
#' @format A `data.table` with "sigminer.features" class name
#' @source Generate from code under data_raw/
#' @examples
#' data(CN.features)
NULL

#' A List of Simulated SBS-96 Catalog Matrix
#'
#' Data from \doi{10.1038/s43018-020-0027-5}.
#' 5 simulated mutation catalogs are used by the paper but only 4 are available.
#' The data are simulated from COSMIC mutational signatures 1, 2, 3, 5, 6, 8,
#' 12, 13, 17 and 18. Each sample is a linear combination of 5 randomly selected
#' signatures with the addiction of Poisson noise. The number of mutation in
#' each sample is randomly selected between 1,000 and 50,000 mutations, in log
#' scale so that a lower number of mutations is more likely to be selected.
#' The proportion of each signature in each sample is also random.
#' @docType data
#' @name simulated_catalogs
#' @format A list of matrix
#' @source Generate from code under data_raw/
#' @examples
#' data(simulated_catalogs)
NULL

generate_matrix_SBS <- function(query, ref_genome) {
  ## TODO handle 4 transcriptional bias categories
  query <- query[query$Variant_Type == "SNP"]
  if (nrow(query) == 0) {
    stop("Zero SNPs to analyze!")
  }

  extract.tbl <- data.table::data.table(
    Chromosome = query$Chromosome, Start = query$Start_Position - 2, End = query$End_Position + 2,
    Reference_Allele = query$Reference_Allele, Tumor_Seq_Allele2 = query$Tumor_Seq_Allele2,
    Tumor_Sample_Barcode = query$Tumor_Sample_Barcode, upstream = query$Start_Position - 20,
    downstream = query$End_Position + 20
  )

  message("=> Extracting 5' and 3' adjacent bases")
  ss <- BSgenome::getSeq(
    x = ref_genome,
    names = extract.tbl$Chromosome,
    start = extract.tbl$Start,
    end = extract.tbl$End,
    as.character = TRUE
  )

  message("=> Extracting +/- 20bp around mutated bases for background C>T estimation")
  updwn <- BSgenome::getSeq(
    x = ref_genome, names = extract.tbl$Chromosome, start = extract.tbl$upstream,
    end = extract.tbl$downstream, as.character = FALSE
  )
  updwn.alphFreq <- data.table::as.data.table(BSgenome::alphabetFrequency(x = updwn))[, c("A", "T", "G", "C")] # Nucleotide frequency
  updwn.tnmFreq <- data.table::as.data.table(Biostrings::trinucleotideFrequency(x = updwn, step = 1))

  extract.tbl$pentanucleotide <- as.character(ss)
  extract.tbl$trinucleotide <- substr(extract.tbl$pentanucleotide, 2, 4)
  extract.tbl$updown <- as.character(updwn)

  extract.tbl <- cbind(extract.tbl, updwn.alphFreq)
  extract.tbl <- cbind(extract.tbl, updwn.tnmFreq[, c("TCA", "TCT", "AGA", "TGA")])
  extract.tbl[, tcw := rowSums(extract.tbl[, c("TCA", "TCT")])]
  extract.tbl[, wga := rowSums(extract.tbl[, c("TGA", "AGA")])]

  ## All combinations
  extract.tbl[, Substitution := paste(extract.tbl$Reference_Allele, extract.tbl$Tumor_Seq_Allele2, sep = ">")]
  extract.tbl$SubstitutionMotif <- paste0(
    substr(x = as.character(extract.tbl$pentanucleotide), 1, 2),
    "[", extract.tbl$Substitution, "]",
    substr(as.character(extract.tbl$pentanucleotide), 4, 5)
  )
  extract.tbl$TriSubstitutionMotif <- substr(extract.tbl$SubstitutionMotif, 2, 8)

  # substitutions are referred to by the pyrimidine of the mutated Watson-Crick base pair
  conv <- c("T>C", "T>C", "C>T", "C>T", "T>A", "T>A", "T>G", "T>G", "C>A", "C>A", "C>G", "C>G")
  names(conv) <- c("A>G", "T>C", "C>T", "G>A", "A>T", "T>A", "A>C", "T>G", "C>A", "G>T", "C>G", "G>C")
  complement <- c("A", "C", "G", "T")
  names(complement) <- c("T", "G", "C", "A")


  extract.tbl$SubstitutionType <- conv[extract.tbl$Substitution]
  # need to reverse-complement triplet for mutated purines (not just the middle base)
  extract.tbl$should_reverse <- extract.tbl$Substitution != extract.tbl$SubstitutionType
  extract.tbl$SubstitutionTypeMotif <- ifelse(extract.tbl$should_reverse,
    paste0(
      complement[substr(x = extract.tbl$pentanucleotide, 5, 5)],
      complement[substr(x = extract.tbl$pentanucleotide, 4, 4)],
      "[", extract.tbl$SubstitutionType, "]",
      complement[substr(x = extract.tbl$pentanucleotide, 2, 2)],
      complement[substr(x = extract.tbl$pentanucleotide, 1, 1)]
    ),
    paste0(
      substr(x = as.character(extract.tbl$pentanucleotide), 1, 2),
      "[", extract.tbl$SubstitutionType, "]",
      substr(as.character(extract.tbl$pentanucleotide), 4, 5)
    )
  )
  extract.tbl$TriSubstitutionTypeMotif <- substr(extract.tbl$SubstitutionTypeMotif, 2, 8)


  # Possible substitution types after being referred to by the pyrimidine of the mutated Watson-Crick base pair
  penta_comb <- expand.grid(
    complement,
    complement,
    "[",
    unique(as.character(conv)),
    "]",
    complement,
    complement,
    stringsAsFactors = FALSE
  ) %>%
    apply(1, paste0, collapse = "") %>%
    unique()

  tri_comb <- substr(penta_comb, 2, 8) %>%
    unique()

  tri_comb2 <- expand.grid(
    complement,
    "[",
    unique(c(as.character(conv), names(conv))),
    "]",
    complement,
    stringsAsFactors = FALSE
  ) %>%
    apply(1, paste0, collapse = "") %>%
    unique()

  # Set levels for type (mainly component)
  extract.tbl$SubstitutionType <- factor(extract.tbl$SubstitutionType, levels = unique(as.character(conv)))
  extract.tbl$TriSubstitutionTypeMotif <- factor(extract.tbl$TriSubstitutionTypeMotif, levels = tri_comb)
  extract.tbl$SubstitutionTypeMotif <- factor(extract.tbl$SubstitutionTypeMotif, levels = penta_comb)

  extract.tbl$TriSubstitutionMotif <- factor(extract.tbl$TriSubstitutionMotif, levels = tri_comb2)

  # Compile data
  ## This is nucleotide frequcny and motif frequency across 41 bp in C>T and C>G context.
  apobecSummary <- extract.tbl[
    as.character(extract.tbl$SubstitutionType) %in% c("C>T", "C>G"),
    .(
      A = sum(A), T = sum(T), G = sum(G), C = sum(C), tcw = sum(tcw),
      wga = sum(wga), bases = sum(A, T, G, C)
    ), Tumor_Sample_Barcode
  ]

  ## This is per sample conversion events
  sub.tbl <- extract.tbl[, .N, list(Tumor_Sample_Barcode, Substitution)]
  sub.tbl <- data.table::dcast(data = sub.tbl, formula = Tumor_Sample_Barcode ~ Substitution, fill = 0, value.var = "N", drop = FALSE)
  sub.tbl[, n_A := rowSums(sub.tbl[, c("A>C", "A>G", "A>T")], na.rm = TRUE)][, n_T := rowSums(sub.tbl[, c("T>A", "T>C", "T>G")], na.rm = TRUE)][, n_G := rowSums(sub.tbl[, c("G>A", "G>C", "G>T")], na.rm = TRUE)][, n_C := rowSums(sub.tbl[, c("C>A", "C>G", "C>T")], na.rm = TRUE)]
  sub.tbl[, n_mutations := rowSums(sub.tbl[, c("n_A", "n_T", "n_G", "n_C")], na.rm = TRUE)]
  sub.tbl[, "n_C>G_and_C>T" := rowSums(sub.tbl[, c("C>G", "G>C", "C>T", "G>A")], na.rm = TRUE)] # number of APOBEC type mutations (C>G and C>T type)

  ## This is per substitution type events
  subType.tbl <- extract.tbl[, .N, .(Tumor_Sample_Barcode, TriSubstitutionMotif)]
  subType.tbl <- data.table::dcast(data = subType.tbl, formula = Tumor_Sample_Barcode ~ TriSubstitutionMotif, fill = 0, value.var = "N", drop = FALSE)

  ### tCw events
  subType.tbl[, tCw_to_A := rowSums(subType.tbl[, .(`T[C>A]A`, `T[C>A]T`)], na.rm = TRUE)]
  subType.tbl[, tCw_to_G := rowSums(subType.tbl[, .(`T[C>G]A`, `T[C>G]T`)], na.rm = TRUE)]
  subType.tbl[, tCw_to_T := rowSums(subType.tbl[, .(`T[C>T]A`, `T[C>T]T`)], na.rm = TRUE)]
  subType.tbl[, tCw := rowSums(subType.tbl[, .(tCw_to_A, tCw_to_G, tCw_to_T)], na.rm = TRUE)]

  ### wGa events
  subType.tbl[, wGa_to_C := rowSums(subType.tbl[, .(`A[G>C]A`, `T[G>C]A`)], na.rm = TRUE)]
  subType.tbl[, wGa_to_T := rowSums(subType.tbl[, .(`A[G>T]A`, `T[G>T]A`)], na.rm = TRUE)]
  subType.tbl[, wGa_to_A := rowSums(subType.tbl[, .(`A[G>A]A`, `T[G>A]A`)], na.rm = TRUE)]
  subType.tbl[, wGa := rowSums(subType.tbl[, .(wGa_to_C, wGa_to_T, wGa_to_A)], na.rm = TRUE)]

  ## tCw_to_G+tCw_to_T
  subType.tbl[, "tCw_to_G+tCw_to_T" := rowSums(subType.tbl[, .(`T[C>G]T`, `T[C>G]A`, `T[C>T]T`, `T[C>T]A`, `T[G>C]A`, `A[G>C]A`, `T[G>A]A`, `A[G>A]A`)], na.rm = TRUE)]

  ### Merge data
  sub.tbl <- merge(sub.tbl, subType.tbl[, .(
    tCw_to_A, tCw_to_T, tCw_to_G, tCw, wGa_to_C, wGa_to_T, wGa_to_A,
    wGa, `tCw_to_G+tCw_to_T`, Tumor_Sample_Barcode
  )],
  by = "Tumor_Sample_Barcode"
  )
  sub.tbl <- merge(sub.tbl, apobecSummary, by = "Tumor_Sample_Barcode")

  ### Estimate APOBEC enrichment
  sub.tbl[, APOBEC_Enrichment := (`tCw_to_G+tCw_to_T` / `n_C>G_and_C>T`) / ((tcw + wga) / (C + G))]
  sub.tbl[, non_APOBEC_mutations := n_mutations - `tCw_to_G+tCw_to_T`]
  sub.tbl[, fraction_APOBEC_mutations := round((n_mutations - non_APOBEC_mutations) / n_mutations, digits = 3)]
  data.table::setDF(sub.tbl)

  message("=> Estimating APOBEC enrichment scores")
  apobec.fisher.dat <- sub.tbl[, c(19, 28, 32, 33, 34)]
  if (nrow(apobec.fisher.dat) == 1) {
    apobec.fisher.dat <- t(as.matrix(apply(X = apobec.fisher.dat, 2, as.numeric)))
  } else {
    apobec.fisher.dat <- apply(X = apobec.fisher.dat, 2, as.numeric)
  }

  ### One way Fisher test to estimate over representation og APOBEC associated tcw mutations
  message("=> Performing one-way Fisher's test for APOBEC enrichment")
  sub.tbl <- cbind(sub.tbl, data.table::rbindlist(apply(X = apobec.fisher.dat, 1, function(x) {
    xf <- fisher.test(matrix(c(x[2], sum(x[3], x[4]), x[1] - x[2], x[3] - x[4]), nrow = 2), alternative = "g")
    data.table::data.table(fisher_pvalue = xf$p.value, or = xf$estimate, ci.up = xf$conf.int[1], ci.low = xf$conf.int[2])
  })))

  data.table::setDT(sub.tbl)
  colnames(sub.tbl)[29:35] <- paste0("n_bg_", colnames(sub.tbl)[29:35])
  sub.tbl <- sub.tbl[order(sub.tbl$fisher_pvalue)]

  ## Choosing APOBEC Enrichment scores > 2 as cutoff
  sub.tbl$APOBEC_Enriched <- ifelse(test = sub.tbl$APOBEC_Enrichment > 2, yes = "yes", no = "no")
  sub.tbl[, fdr := p.adjust(sub.tbl$fisher_pvalue, method = "fdr")] # Adjusted p-values

  message(
    paste0("=> APOBEC related mutations are enriched in "),
    round(nrow(sub.tbl[APOBEC_Enriched %in% "yes"]) / nrow(sub.tbl) * 100, digits = 3),
    "% of samples (APOBEC enrichment score > 2; ",
    nrow(sub.tbl[APOBEC_Enriched %in% "yes"]), " of ", nrow(sub.tbl), " samples)"
  )

  message("=> Creating SBS sample-by-component matrice")

  SBS_6 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "SubstitutionType")
  SBS_6 <- SBS_6[, c("T>C", "C>T", "T>A", "T>G", "C>A", "C>G")] %>% as.matrix()
  SBS_96 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "TriSubstitutionTypeMotif")
  SBS_96 <- SBS_96[, tri_comb] %>% as.matrix()
  SBS_1536 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "SubstitutionTypeMotif")
  SBS_1536 <- SBS_1536[, penta_comb] %>% as.matrix()

  res <- list(
    nmf_matrix = SBS_96,
    all_matrice = list(
      SBS_6 = SBS_6,
      SBS_96 = SBS_96,
      SBS_1536 = SBS_1536
    ),
    APOBEC_scores = sub.tbl
  )
  res
}

generate_matrix_INDEL <- function(query, ref_genome, genome_build = "hg19", add_trans_bias = FALSE) {
  query <- subset(query,query$Variant_Type != "SNP")
  if (nrow(query) == 0) {
    stop("Zero INDELs to analyze!")
  }
  indel_types = c("1:Del:C:0", "1:Del:C:1", "1:Del:C:2", "1:Del:C:3", "1:Del:C:4", "1:Del:C:5",
                  "1:Del:T:0", "1:Del:T:1", "1:Del:T:2", "1:Del:T:3", "1:Del:T:4", "1:Del:T:5",
                  "1:Ins:C:0", "1:Ins:C:1", "1:Ins:C:2", "1:Ins:C:3", "1:Ins:C:4", "1:Ins:C:5",
                  "1:Ins:T:0", "1:Ins:T:1", "1:Ins:T:2", "1:Ins:T:3", "1:Ins:T:4", "1:Ins:T:5",
                  # >1bp INDELS
                  "2:Del:R:0", "2:Del:R:1", "2:Del:R:2", "2:Del:R:3", "2:Del:R:4", "2:Del:R:5",
                  "3:Del:R:0", "3:Del:R:1", "3:Del:R:2", "3:Del:R:3", "3:Del:R:4", "3:Del:R:5",
                  "4:Del:R:0", "4:Del:R:1", "4:Del:R:2", "4:Del:R:3", "4:Del:R:4", "4:Del:R:5",
                  "5:Del:R:0", "5:Del:R:1", "5:Del:R:2", "5:Del:R:3", "5:Del:R:4", "5:Del:R:5",
                  "2:Ins:R:0", "2:Ins:R:1", "2:Ins:R:2", "2:Ins:R:3", "2:Ins:R:4", "2:Ins:R:5",
                  "3:Ins:R:0", "3:Ins:R:1", "3:Ins:R:2", "3:Ins:R:3", "3:Ins:R:4", "3:Ins:R:5",
                  "4:Ins:R:0", "4:Ins:R:1", "4:Ins:R:2", "4:Ins:R:3", "4:Ins:R:4", "4:Ins:R:5",
                  "5:Ins:R:0", "5:Ins:R:1", "5:Ins:R:2", "5:Ins:R:3", "5:Ins:R:4", "5:Ins:R:5",
                  #MicroHomology INDELS
                  "2:Del:M:1", "3:Del:M:1", "3:Del:M:2", "4:Del:M:1", "4:Del:M:2", "4:Del:M:3",
                  "5:Del:M:1", "5:Del:M:2", "5:Del:M:3", "5:Del:M:4", "5:Del:M:5",
                  "complex", "non_matching")
  # 初始化矩阵和内容
  indel_type_count <- 0
  rowname <- indel_types
  colname <- query$Tumor_Sample_Barcode
  ID_matrix <- matrix(data = 0,nrow = 85, ncol = nrow(query), byrow = FALSE, dimnames = list(rowname,colname))
  for (i in 1:nrow(query)){
    print(i)
    data_sample <- query[c(i),]
    # indel_type_count <- paste(indel_type_count,0)
    ref <- data_sample$Reference_Allele
    mut <- data_sample$Tumor_Seq_Allele2
    mutation_ID <- data.table(indel_types,indel_type_count)
    indel_type_count <- 0
    ###### Tips:为了减少条件语句的使用，将该样本突变处有序列的部分赋值给ID_type
    if (ref != '-' & mut == '-'){
      ID_type <- ref
    }else if(mut != '-' & ref == '-'){
      ID_type <- mut
    }
    ###### 先将complex分出来
    if(ref != '-' & mut != '-'){
      real_type <- 'complex'
      for (i in 1:96){
        if(real_type == indel_types[i]){
          mutation_ID$indel_type_count[i] <- mutation_ID$indel_type_count[i] + 1
          break
        }
      }
    }else{
      ###### Part2.1 记录mutation类型是ins还是del
      if (mut == '-'){
        mut_type = 'Del'
        mut_base = substr(ref,1,1)
      }else{
        mut_type = 'Ins'
        mut_base = substr(mut,1,1)
      }
      ###### Part3. 分类并计数indel_type
      ###### Part3.1 取上下游序列
      Chromosome=paste0("chr",data_sample$Chromosome)
      Start=data_sample$Start_Position-50
      End = data_sample$End_Position+250
      data_loc <- data.frame(Chromosome,Start,End)
      # ref_genome <- BSgenome.Hsapiens.UCSC.hg19
      get_seq <- BSgenome::getSeq(x = ref_genome,
                                  names = data_loc$Chromosome,
                                  start = data_loc$Start,
                                  end = data_loc$End)
      upstream = BSgenome::getSeq(x = ref_genome,
                                  names = data_loc$Chromosome,
                                  start = data_loc$Start,
                                  end = data_loc$Start+49)
      downstream = BSgenome::getSeq(x = ref_genome,
                                    names = data_loc$Chromosome,
                                    start = data_loc$End-249,
                                    end = data_loc$End)
      ###### Part3.2 计数突变处序列长度:length_INDEL
      length_INDEL <- 0
      if (nchar(as.character(ID_type)) < 5){
        length_INDEL <- nchar(as.character(ID_type))
      }else if(nchar(as.character(ID_type)) >= 5){
        length_INDEL <- 5
      }
      ###### Part3.3 记录和突变位点匹配的次数
      # a = "AA"
      # # b = "AAAACCTT"
      # paste0()
      # sub("(^AA+).*", "\\1", b)
      # sub("^((AA)+).*$", "\\1", "AAAACC")
      count_repeat = function(x, sequence) {
        if (startsWith(sequence, x)) {
          pattern = paste0(
            "^((",
            x,
            ")+).*$"
          )
          nchar(sub(pattern,"\\1", sequence)) / nchar(x)
        } else {
          return(0L)
        }
      }
      count_repeat <- count_repeat(as.character(ID_type),as.character(downstream))
      ###### Part3.4 记录homology次数
      count_homology_size = function(x, upstream, downstream) {
        size = nchar(as.character(x))
        if (size >= 2) {
          x_down = substr(x, 1, size - 1)
          size_down = 0
          while(nchar(x_down) > 0) {
            if (startsWith(downstream, x_down)) {
              size_down = nchar(x_down)
              break()
            } else {
              x_down = substr(x_down, 1, nchar(x_down) -1)
            }
          }
          # Upstream
          # ACAAC|TC|AAGCGGC
          x_up = substring(x, 2)
          size_up = 0
          while(nchar(x_up) > 0) {
            if (endsWith(upstream, x_up)) {
              size_up = nchar(x_up)
              break()
            } else {
              x_up = substring(x_up, 2)
            }
          }
          max(size_up, size_down)
        }else{
          0
        }
      }
      count_homology_size <- count_homology_size(ID_type,as.character(upstream), as.character(downstream))
      ###### Part3.4 判断单位点还是多位点，区别计数
      ###### repeat_num是重复数的计数（单个位点和多个位点的homology计数方法不同）
      repeat_num <- 0
      ###### Part3.4.1 单个位点
      # 转化整个链（针对1bp indel）
      if (nchar(as.character(ID_type)) == 1){
        if (mut_base == 'G'){
          mut_base = 'C'
          strand <- '-1'
        }else if(mut_base == 'A'){
          mut_base = 'T'
          strand <- '-1'
        }
        if(count_repeat > 5){
          count_repeat <- 5
        }
        ###### Part3.4.2 多个位点
      }else if(nchar(as.character(ID_type)) > 1){
        mut_base <- 'R'
        if(count_repeat > 5){
          count_repeat <- 5
        }
      }else if(count_homology_size > 0){
        mut_base <- 'M'
        if(count_homology_size > 5){
          count_homology_size <- 5
        }
      }
      ###### Part4 整合indel_type为标准格式，记录为real_types
      if(nchar(ID_type == 1)){
        real_type <- paste(length_INDEL,":",mut_type,":",mut_base,":",count_repeat,sep="")
      }else if(mut_base == 'R' & mut_type == 'Del'){
        real_type <- paste(length_INDEL,":",mut_type,":",mut_base,":",count_repeat+1,sep="")
      }else if(mut_base == 'R' & mut_type == 'Ins'){
        real_type <- paste(length_INDEL,":",mut_type,":",mut_base,":",count_repeat,sep="")
      }else{
        real_type <- paste(length_INDEL,":",mut_type,":",mut_base,":",count_homology_size,sep="")
      }
      ###### Part5 对indel_type进行计数
      indel_type_count <- 0
      mutation_ID <- data.table(indel_types,indel_type_count)
      for (n in 1:85){
        if(identical(real_type,indel_types[n])){
          # mutation_ID$indel_type_count[n] <- mutation_ID$indel_type_count[n] + 1
          ID_matrix[n,i] <- ID_matrix[n,i] + 1
        }
      }
    }
    # 检查单个样本的类型：mutation_ID[indel_type_count > 0]
  }
}

records_to_matrix <- function(dt, samp_col, component_col) {
  dt.summary <- dt[, .N, by = c(samp_col, component_col)]
  mat <- as.data.frame(data.table::dcast(dt.summary,
    formula = as.formula(paste(samp_col, "~", component_col)),
    fill = 0, value.var = "N", drop = FALSE
  ))
  rownames(mat) <- mat[, 1]
  mat <- mat[, -1]
  mat
}

utils::globalVariables(
  c("A", "APOBEC_Enriched", "APOBEC_Enrichment",
    "A[G>A]A", "A[G>C]A", "A[G>T]A", "G",
    "Substitution",
    "T[C>A]A", "T[C>A]T",
    "T[C>G]A", "T[C>G]T",
    "T[C>T]A", "T[C>T]T", "T[G>A]A",
    "T[G>C]A", "T[G>T]A",
    "TriSubstitutionMotif", "Tumor_Sample_Barcode", "fdr",
    "fraction_APOBEC_mutations", "n_A", "n_C", "n_C>G_and_C>T",
    "n_G", "n_T", "n_mutations", "non_APOBEC_mutations",
    "tCw", "tCw_to_A", "tCw_to_G", "tCw_to_G+tCw_to_T",
    "tCw_to_T", "tcw", "wGa", "wGa_to_A", "wGa_to_C","wGa_to_T", "wga")
)

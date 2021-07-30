
# SBS ---------------------------------------------------------------------

generate_matrix_SBS <- function(query, ref_genome, genome_build = "hg19", add_trans_bias = FALSE) {
  send_info("SBS matrix generation - start.")
  ## T: and U: are basically same to sigprofiler, the other two types (N: and B:) are
  ## a bit different, may be caused by different annotations.
  query <- query[query$Variant_Type == "SNP"]
  if (nrow(query) == 0) {
    stop("Zero SNPs to analyze!")
  }

  extract.tbl <- data.table::data.table(
    Chromosome = query$Chromosome, Start = query$Start_Position, End = query$End_Position,
    Reference_Allele = query$Reference_Allele, Tumor_Seq_Allele2 = query$Tumor_Seq_Allele2,
    Tumor_Sample_Barcode = query$Tumor_Sample_Barcode, upstream = query$Start_Position - 20,
    downstream = query$End_Position + 20
  )

  send_info("Extracting 5' and 3' adjacent bases.")
  ss <- BSgenome::getSeq(
    x = ref_genome,
    names = extract.tbl$Chromosome,
    start = extract.tbl$Start - 2,
    end = extract.tbl$End + 2,
    as.character = TRUE
  )

  send_info("Extracting +/- 20bp around mutated bases for background C>T estimation.")
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

  # penta_comb <- expand.grid(
  #   complement,
  #   complement,
  #   "[",
  #   unique(as.character(conv)),
  #   "]",
  #   complement,
  #   complement,
  #   stringsAsFactors = FALSE
  # ) %>%
  #   apply(1, paste0, collapse = "") %>%
  #   unique()

  penta_comb <- vector_to_combination(
    complement,
    complement,
    "[",
    unique(as.character(conv)),
    "]",
    complement,
    complement
  )

  tri_comb <- substr(penta_comb, 2, 8) %>%
    unique()

  tri_comb2 <- vector_to_combination(
    complement,
    "[",
    unique(c(as.character(conv), names(conv))),
    "]",
    complement
  )

  # Set levels for type (mainly component)
  extract.tbl$SubstitutionType <- factor(extract.tbl$SubstitutionType, levels = unique(as.character(conv)))
  extract.tbl$TriSubstitutionTypeMotif <- factor(extract.tbl$TriSubstitutionTypeMotif, levels = tri_comb)
  extract.tbl$SubstitutionTypeMotif <- factor(extract.tbl$SubstitutionTypeMotif, levels = penta_comb)

  extract.tbl$TriSubstitutionMotif <- factor(extract.tbl$TriSubstitutionMotif, levels = tri_comb2)
  extract.tbl$Substitution <- factor(extract.tbl$Substitution, levels = unique(c(as.character(conv), names(conv))))

  # Compile data
  ## This is nucleotide frequcny and motif frequency across 41 bp in C>T and C>G context.
  apobecSummary <- extract.tbl[
    as.character(extract.tbl$SubstitutionType) %in% c("C>T", "C>G"),
    .(
      A = sum(A), T = sum(T), G = sum(G), C = sum(C), tcw = sum(tcw),
      wga = sum(wga), bases = sum(A, T, G, C)
    ), Tumor_Sample_Barcode
  ]
  ## The by operation may remove some samples here
  ## should it be added??

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

  send_info("Estimating APOBEC enrichment scores.")
  apobec.fisher.dat <- sub.tbl[, c(19, 28, 32, 33, 34)]
  if (nrow(apobec.fisher.dat) == 1) {
    apobec.fisher.dat <- t(as.matrix(apply(X = apobec.fisher.dat, 2, as.numeric)))
  } else {
    apobec.fisher.dat <- apply(X = apobec.fisher.dat, 2, as.numeric)
  }

  ### One way Fisher test to estimate over representation og APOBEC associated tcw mutations
  send_info("Performing one-way Fisher's test for APOBEC enrichment.")
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

  send_success(
    paste0("APOBEC related mutations are enriched in "),
    round(nrow(sub.tbl[APOBEC_Enriched %in% "yes"]) / nrow(sub.tbl) * 100, digits = 3),
    "% of samples (APOBEC enrichment score > 2; ",
    nrow(sub.tbl[APOBEC_Enriched %in% "yes"]), " of ", nrow(sub.tbl), " samples)"
  )

  send_info("Creating SBS sample-by-component matrices.")

  SBS_6 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "SubstitutionType")
  SBS_6 <- SBS_6[, c("T>C", "C>T", "T>A", "T>G", "C>A", "C>G")] %>% as.matrix()
  send_success("SBS-6 matrix created.")

  SBS_96 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "TriSubstitutionTypeMotif")
  SBS_96 <- SBS_96[, tri_comb] %>% as.matrix()
  send_success("SBS-96 matrix created.")

  SBS_1536 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "SubstitutionTypeMotif")
  SBS_1536 <- SBS_1536[, penta_comb] %>% as.matrix()
  send_success("SBS-1536 matrix created.")

  if (add_trans_bias) {
    t_labels <- c("T:", "U:", "B:", "N:")

    SBS_24 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "SubstitutionType",
      add_trans_bias = TRUE, build = genome_build
    )
    SBS_24 <- SBS_24[, vector_to_combination(
      t_labels,
      c("T>C", "C>T", "T>A", "T>G", "C>A", "C>G")
    )] %>% as.matrix()
    send_success("SBS-24 (6x4) matrix created.")

    SBS_384 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "TriSubstitutionTypeMotif",
      add_trans_bias = TRUE, build = genome_build
    )
    SBS_384 <- SBS_384[, vector_to_combination(t_labels, tri_comb)] %>% as.matrix()
    send_success("SBS-384 (96x4) matrix created.")

    SBS_6144 <- records_to_matrix(extract.tbl, "Tumor_Sample_Barcode", "SubstitutionTypeMotif",
      add_trans_bias = TRUE, build = genome_build
    )
    SBS_6144 <- SBS_6144[, vector_to_combination(t_labels, penta_comb)] %>% as.matrix()
    send_success("SBS-6144 (1536x4) matrix created.")
  }

  if (add_trans_bias) {
    send_info("Return SBS-192 as major matrix.")
    SBS_192 <- SBS_384[, grepl("T:|U:", colnames(SBS_384))]
    res <- list(
      nmf_matrix = SBS_192,
      all_matrices = list(
        SBS_6 = SBS_6,
        SBS_24 = SBS_24,
        SBS_96 = SBS_96,
        SBS_192 = SBS_192,
        SBS_384 = SBS_384,
        SBS_1536 = SBS_1536,
        SBS_6144 = SBS_6144
      ),
      APOBEC_scores = sub.tbl
    )
  } else {
    send_info("Return SBS-96 as major matrix.")
    res <- list(
      nmf_matrix = SBS_96,
      all_matrices = list(
        SBS_6 = SBS_6,
        SBS_96 = SBS_96,
        SBS_1536 = SBS_1536
      ),
      APOBEC_scores = sub.tbl
    )
  }
  res
}


# DBS ---------------------------------------------------------------------

generate_matrix_DBS <- function(query, ref_genome, genome_build = "hg19", add_trans_bias = FALSE) {
  send_info("DBS matrix generation - start.")

  query <- query[query$Variant_Type %in% c("SNP", "DNP")]
  if (nrow(query) == 0) {
    send_stop("Zero SNP/DNPs to analyze!")
  }

  ## Generate DBS catalogues
  conv <- c(
    "AC>CA", "AC>CG", "AC>CT", "AC>GA", "AC>GG", "AC>GT",
    "AC>TA", "AC>TG", "AC>TT", "AT>CA", "AT>CC", "AT>CG",
    "AT>GA", "AT>GC", "AT>TA", "CC>AA", "CC>AG", "CC>AT",
    "CC>GA", "CC>GG", "CC>GT", "CC>TA", "CC>TG", "CC>TT",
    "CG>AT", "CG>GC", "CG>GT", "CG>TA", "CG>TC", "CG>TT",
    "CT>AA", "CT>AC", "CT>AG", "CT>GA", "CT>GC", "CT>GG",
    "CT>TA", "CT>TC", "CT>TG", "GC>AA", "GC>AG", "GC>AT",
    "GC>CA", "GC>CG", "GC>TA", "TA>AT", "TA>CG", "TA>CT",
    "TA>GC", "TA>GG", "TA>GT", "TC>AA", "TC>AG", "TC>AT",
    "TC>CA", "TC>CG", "TC>CT", "TC>GA", "TC>GG", "TC>GT",
    "TG>AA", "TG>AC", "TG>AT", "TG>CA", "TG>CC", "TG>CT",
    "TG>GA", "TG>GC", "TG>GT", "TT>AA", "TT>AC", "TT>AG",
    "TT>CA", "TT>CC", "TT>CG", "TT>GA", "TT>GC", "TT>GG"
  )
  names(conv) <- c(
    "GT>TG", "GT>CG", "GT>AG", "GT>TC", "GT>CC", "GT>AC",
    "GT>TA", "GT>CA", "GT>AA", "AT>TG", "AT>GG", "AT>CG",
    "AT>TC", "AT>GC", "AT>TA", "GG>TT", "GG>CT", "GG>AT",
    "GG>TC", "GG>CC", "GG>AC", "GG>TA", "GG>CA", "GG>AA",
    "CG>AT", "CG>GC", "CG>AC", "CG>TA", "CG>GA", "CG>AA",
    "AG>TT", "AG>GT", "AG>CT", "AG>TC", "AG>GC", "AG>CC",
    "AG>TA", "AG>GA", "AG>CA", "GC>TT", "GC>CT", "GC>AT",
    "GC>TG", "GC>CG", "GC>TA", "TA>AT", "TA>CG", "TA>AG",
    "TA>GC", "TA>CC", "TA>AC", "GA>TT", "GA>CT", "GA>AT",
    "GA>TG", "GA>CG", "GA>AG", "GA>TC", "GA>CC", "GA>AC",
    "CA>TT", "CA>GT", "CA>AT", "CA>TG", "CA>GG", "CA>AG",
    "CA>TC", "CA>GC", "CA>AC", "AA>TT", "AA>GT", "AA>CT",
    "AA>TG", "AA>GG", "AA>CG", "AA>TC", "AA>GC", "AA>CC"
  )

  ## Search for DBS
  send_info("Searching DBS records...")
  query_DNP <- query[
    query$Variant_Type == "DNP",
    c(
      "Tumor_Sample_Barcode", "Hugo_Symbol", "Chromosome",
      "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"
    )
  ]
  query <- query[, search_DBS(.SD), by = Tumor_Sample_Barcode]
  query <- rbind(query, query_DNP, fill = TRUE)
  send_success("Done.")

  if (nrow(query) == 0) {
    send_stop("Zero DBSs to analyze!")
  }

  query[, dbs := paste(Reference_Allele, Tumor_Seq_Allele2, sep = ">")]
  query[, dbsMotif := ifelse(dbs %in% names(conv), conv[dbs], dbs)]
  query[, dbsMotif := factor(dbsMotif, levels = as.character(conv))]
  query[, should_reverse := substr(dbsMotif, 1, 2) != substr(dbs, 1, 2)]

  ## Query sequence
  query[, upstream := as.character(
    BSgenome::getSeq(
      x = ref_genome,
      names = Chromosome,
      start = Start_Position - 1,
      end = Start_Position - 1
    )
  )]
  query[, downstream := as.character(
    BSgenome::getSeq(
      x = ref_genome,
      names = Chromosome,
      start = Start_Position + 2,
      end = Start_Position + 2
    )
  )]
  send_success("Reference sequences queried from genome.")

  query$complexMotif <- paste0(query$upstream, "[", query$dbsMotif, "]", query$downstream)
  query$complexMotif <- factor(query$complexMotif,
    levels = vector_to_combination(
      c("A", "C", "G", "T"),
      "[",
      levels(query$dbsMotif),
      "]",
      c("A", "C", "G", "T")
    )
  )

  DBS_78 <- records_to_matrix(query, "Tumor_Sample_Barcode", "dbsMotif") %>% as.matrix()
  send_success("DBS-78 matrix created.")

  DBS_1248 <- records_to_matrix(query, "Tumor_Sample_Barcode", "complexMotif") %>% as.matrix()
  send_success("DBS-1248 matrix created.")

  if (add_trans_bias) {
    query$End_Position <- query$Start_Position + 1
    DBS_186 <- records_to_matrix(query, "Tumor_Sample_Barcode", "dbsMotif",
      add_trans_bias = TRUE, build = genome_build, mode = "DBS"
    )
    DBS_186 <- DBS_186 %>% as.matrix()
    send_success("DBS-186 matrix created.")

    send_info("Return DBS-186 as major matrix.")
    res <- list(
      nmf_matrix = DBS_186,
      all_matrices = list(
        DBS_78 = DBS_78,
        DBS_186 = DBS_186,
        DBS_1248 = DBS_1248
      )
    )
  } else {
    send_info("Return SBS-78 as major matrix.")
    res <- list(
      nmf_matrix = DBS_78,
      all_matrices = list(
        DBS_78 = DBS_78,
        DBS_1248 = DBS_1248
      )
    )
  }
  res
}

## Search DBS in SNV records
search_DBS <- function(x) {
  # x = data.table::as.data.table(x)
  ## Set a position data.table to cover all combinations
  pos_dt <- rbind(
    data.table::data.table(
      chr = x$Chromosome,
      start = x$Start_Position - 1,
      end = x$Start_Position
    ),
    data.table::data.table(
      chr = x$Chromosome,
      start = x$Start_Position,
      end = x$Start_Position
    )
  ) %>%
    unique()

  data.table::setkey(pos_dt, chr, start, end)

  ## Find the overlaps
  x$chr <- x$Chromosome
  x$start <- x$Start_Position
  x$end <- x$start

  overlap_dt <- data.table::foverlaps(x, pos_dt, type = "within")
  overlap_dt[, region := paste0(chr, ":", start, "-", end)]

  freq_dt <- overlap_dt[, .N, by = region]

  DBS_index <- freq_dt$N > 1
  if (sum(DBS_index) < 1) {
    return(data.table::data.table())
  } else {
    DBS_regions <- freq_dt$region[DBS_index]
    res_dt <- overlap_dt[region %in% DBS_regions][order(Start_Position)] # Mutation position should be ordered
    return(
      res_dt[, list(
        Hugo_Symbol = paste(unique(Hugo_Symbol), collapse = ","),
        Chromosome = paste(unique(Chromosome), collapse = ","),
        Start_Position = Start_Position[1],
        Reference_Allele = paste(Reference_Allele, collapse = ""),
        Tumor_Seq_Allele2 = paste(Tumor_Seq_Allele2, collapse = "")
      ), by = region]
    )
  }
}

# INDELs (ID) -------------------------------------------------------------

adjust_indels <- function(query) {
  query %>%
    dplyr::rowwise() %>%
    dplyr::mutate(cN = nchar(Biobase::lcPrefixC(c(.data$Reference_Allele, .data$Tumor_Seq_Allele2), ignore.case = TRUE))) %>%
    dplyr::mutate(
      Start_Position = .data$Start_Position + .data$cN,
      Reference_Allele = substring(.data$Reference_Allele, .data$cN + 1L),
      Tumor_Seq_Allele2 = substring(.data$Tumor_Seq_Allele2, .data$cN + 1L)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      Reference_Allele = ifelse(.data$Reference_Allele == "", "-", .data$Reference_Allele),
      Tumor_Seq_Allele2 = ifelse(.data$Tumor_Seq_Allele2 == "", "-", .data$Tumor_Seq_Allele2),
      End_Position = .data$Start_Position + nchar(.data$Reference_Allele) - 1L
    ) %>%
    data.table::as.data.table()
}

generate_matrix_INDEL <- function(query, ref_genome, genome_build = "hg19", add_trans_bias = FALSE) {
  send_info("INDEL matrix generation - start.")

  query <- query[query$Variant_Type %in% c("INS", "DEL")]
  if (nrow(query) == 0) {
    send_stop("Zero INDELs to analyze!")
  }

  indel_types <- c(
    "1:Del:C:0", "1:Del:C:1", "1:Del:C:2", "1:Del:C:3", "1:Del:C:4", "1:Del:C:5",
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
    # MicroHomology INDELS
    "2:Del:M:1", "3:Del:M:1", "3:Del:M:2", "4:Del:M:1", "4:Del:M:2", "4:Del:M:3",
    "5:Del:M:1", "5:Del:M:2", "5:Del:M:3", "5:Del:M:4", "5:Del:M:5",
    "complex"
    # , "non_matching"
  )

  ## Adjust un-standard variant formats
  query <- adjust_indels(query)
  ## Update Variant_Type
  query[, Variant_Type := ifelse(Variant_Type == "INS", "Ins", "Del")]
  ## Set ID type
  query[, ID_type := ifelse(Reference_Allele == "-",
    Tumor_Seq_Allele2,
    Reference_Allele
  )]
  ## Seach 'complex' motif
  query[, ID_motif := ifelse(Reference_Allele != "-" & Tumor_Seq_Allele2 != "-",
    "complex", NA_character_
  )]

  ## Drop records with 'complex' motif
  query_comp <- query[ID_motif == "complex"]
  query <- query[is.na(ID_motif)]

  if (nrow(query) == 0) {
    send_stop("Zero INDELs to analyze after dropping 'complex' labeled records!")
  }

  ## Get first INDEL base
  query[, mut_base := substr(ID_type, 1, 1)]

  ## Query sequence
  query[, upstream := as.character(
    BSgenome::getSeq(
      x = ref_genome,
      names = Chromosome,
      start = Start_Position - 50,
      end = Start_Position - 1
    )
  )]

  # NOTE: deletion variants should start from Start_Position
  query[, downstream := as.character(
    BSgenome::getSeq(
      x = ref_genome,
      names = Chromosome,
      start = ifelse(Variant_Type == "Ins", Start_Position, End_Position + 1),
      end = ifelse(Variant_Type == "Ins", Start_Position + 249, End_Position + 250)
    )
  )]
  send_success("Reference sequences queried from genome.")

  ## Length of INDEL
  query[, ID_len := ifelse(nchar(ID_type) < 5, nchar(ID_type), 5)]
  send_success("INDEL length extracted.")

  ## Adjacent repeats (copies) count
  query[, count_repeat := mapply(count_repeat, ID_type, downstream) + mapply(count_repeat, ID_type, upstream, is_downstream = FALSE)]
  send_success("Adjacent copies counted.")
  ## Micro-homology count
  query[, count_homosize := mapply(
    count_homology_size,
    ID_type, upstream, downstream
  ) %>% as.integer()]
  send_success("Microhomology size calculated.")

  ## Set maximum repeat/homology-size count
  query[, count_repeat := ifelse(count_repeat > 5, 5, count_repeat)]
  query[, count_homosize := ifelse(count_homosize > 5, 5, count_homosize)]

  ## For length-1 INDEL
  conv <- c("C", "T")
  names(conv) <- c("G", "A")
  query[, should_reverse := ID_len == 1 & mut_base %in% c("G", "A")]
  query[, mut_base := ifelse(should_reverse, conv[mut_base], mut_base)]
  ## For length-n (n>1) INDEL
  query[, mut_base := ifelse(ID_len > 1, "R", mut_base)]

  ## 这里可能有问题 M 标记
  query[, mut_base := ifelse(Variant_Type == "Del" &
                               ID_len > 1 & count_repeat == 0 &
                               count_homosize > 0, "M", mut_base)]

  ## Generate variables to handle strand bias labeling
  query[, should_reverse := sapply(ID_type, is_all_purine)]
  query[, ID_type := ifelse(should_reverse, purine2pyrimidine(ID_type), ID_type)]
  # determine pyrimidine situaiton
  query$all_pyrimidine <- sapply(query$ID_type, is_all_pyrimidine)

  ## Generate all factors
  sim_types <- c(indel_types[1:24], "long_Del", "long_Ins", "MH", "complex")
  query <- rbind(query, query_comp, fill = TRUE)

  ## Generate Motifs
  ## ifelse() has problems for assign values here
  ## use case_when
  ##
  ## NOTE: according to reference paper says
  ## "Since almost no insertions with microhomologies
  ## were identified across more than 20,000 tumors"
  ## so records like this will classified into 'non_matching' in ID_motif
  ## and 'MH' in ID_motif_sp
  query <- query %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(
      ID_motif = dplyr::case_when(
        ID_len == 1 ~ paste(ID_len, Variant_Type, mut_base, count_repeat, sep = ":"),
        mut_base == "R" ~ paste(ID_len, Variant_Type, mut_base, count_repeat, sep = ":"),
        mut_base == "M" ~ paste(ID_len, Variant_Type, mut_base, count_homosize, sep = ":"),
        ID_motif == "complex" ~ "complex",
        TRUE ~ "non_matching"
      ),
      ID_motif_sp = ID_motif
    ) %>%
    dplyr::mutate(
      ID_motif_sp = dplyr::case_when(
        ID_motif_sp %in% c(indel_types[1:24], "complex") ~ ID_motif_sp,
        ID_motif_sp %in% indel_types[25:48] ~ "long_Del",
        ID_motif_sp %in% indel_types[49:72] ~ "long_Ins",
        TRUE ~ "MH"
      )
    ) %>%
    data.table::as.data.table()

  ## "complex" type will be dropped off from ID-83
  query[, ID_motif := factor(ID_motif, levels = indel_types[-84])]
  query[, ID_motif_sp := factor(ID_motif_sp, levels = sim_types)]
  send_success("INDEL records classified into different components (types).")

  ID_28 <- records_to_matrix(query, "Tumor_Sample_Barcode", "ID_motif_sp") %>% as.matrix()
  send_success("ID-28 matrix created.")
  ID_83 <- records_to_matrix(query, "Tumor_Sample_Barcode", "ID_motif") %>% as.matrix()
  send_success("ID-83 matrix created.")

  if (add_trans_bias) {
    ID_415 <- records_to_matrix(query, "Tumor_Sample_Barcode", "ID_motif",
      add_trans_bias = TRUE, build = genome_build, mode = "ID"
    )
    send_success("ID-415 matrix created.")

    send_info("Return ID-415 as major matrix.")
    res <- list(
      nmf_matrix = ID_415,
      all_matrices = list(
        ID_28 = ID_28,
        ID_83 = ID_83,
        ID_415 = ID_415
      )
    )
  } else {
    send_info("Return ID-83 as major matrix.")
    res <- list(
      nmf_matrix = ID_83,
      all_matrices = list(
        ID_28 = ID_28,
        ID_83 = ID_83
      )
    )
  }
  res
}


# Utils -------------------------------------------------------------------

records_to_matrix <- function(dt, samp_col, component_col, add_trans_bias = FALSE, build = "hg19", mode = "SBS") {
  if (add_trans_bias) {
    transcript_dt <- get_genome_annotation(data_type = "transcript", genome_build = build)
    data.table::setkey(transcript_dt, chrom, start, end)

    if ("Start" %in% colnames(dt)) {
      loc_dt <- dt[, .(Chromosome, Start, End)]
    } else {
      loc_dt <- dt[, .(Chromosome, Start_Position, End_Position)]
    }
    colnames(loc_dt)[1:3] <- c("chrom", "start", "end")
    loc_dt$MutIndex <- 1:nrow(loc_dt)

    m_dt <- data.table::foverlaps(loc_dt, transcript_dt, type = "any")[
      , .(MatchCount = .N, strand = paste0(unique(strand), collapse = "/")),
      by = list(MutIndex)
    ]

    ## Actually, the MatchCount should only be 0, 1, 2
    if (any(m_dt$MatchCount > 2)) {
      send_stop("More than 2 regions counted, please report your data and code to developer!")
    }

    if (mode == "SBS") {
      dt$transcript_bias_label <- ifelse(
        m_dt$MatchCount == 2, "B:",
        ifelse(m_dt$MatchCount == 0, "N:",
          ifelse(xor(dt$should_reverse, m_dt$strand == "-"),
            # (dt$should_reverse & m_dt$strand == "+") | (!dt$should_reverse & m_dt$strand == "-"),
            "T:", "U:"
          )
        ) ## If should reverse, the base switch to template strand from coding strand
      )
      new_levels <- vector_to_combination(
        c("T:", "U:", "B:", "N:"),
        levels(dt[[component_col]])
      )
      dt[[component_col]] <- paste0(dt$transcript_bias_label, dt[[component_col]])
      dt[[component_col]] <- factor(dt[[component_col]], levels = new_levels)
    } else if (mode == "DBS") {

      # dt <- dt %>%
      #   dplyr::as_tibble() %>%
      #   #dplyr::bind_cols(m_dt %>% dplyr::as_tibble()) %>%
      #   dplyr::mutate(
      #     transcript_bias_label = dplyr::case_when(
      #       !substr(.data[[component_col]], 1, 2) %in% c("TT", "TC", "CT", "CC") ~ "Q:",
      #       m_dt$MatchCount == 2 ~ "B:",
      #       m_dt$MatchCount == 0 ~ "N:",
      #       xor(dt$should_reverse, m_dt$strand == "-") ~ "T:",
      #       TRUE ~ "U:"
      #     )
      #   ) %>%
      #   data.table::as.data.table()

      dt$transcript_bias_label <- ifelse(
        substr(dt[[component_col]], 1, 2) %in% c("TT", "TC", "CT", "CC"),
        ifelse(
          m_dt$MatchCount == 2, "B:",
          ifelse(m_dt$MatchCount == 0, "N:",
            ifelse(xor(dt$should_reverse, m_dt$strand == "-"),
              "T:", "U:"
            )
          )
        ), "Q:"
      )

      new_levels <- vector_to_combination(
        c("T:", "U:", "B:", "N:", "Q:"),
        levels(dt[[component_col]])
      )
      new_levels <- ifelse(substr(new_levels, 3, 4) %in% c("TT", "TC", "CT", "CC") & substr(new_levels, 1, 1) != "Q", new_levels,
        ifelse(substr(new_levels, 3, 4) %in% c("TT", "TC", "CT", "CC") & substr(new_levels, 1, 1) == "Q",
          NA, gsub("[A-Z]\\:", "Q:", new_levels)
        )
      ) %>%
        na.omit() %>%
        unique()

      dt[[component_col]] <- paste0(dt$transcript_bias_label, dt[[component_col]])
      dt[[component_col]] <- factor(dt[[component_col]], levels = new_levels)
    } else if (mode == "ID") {
      dt$transcript_bias_label <- ifelse(
        dt$all_pyrimidine == TRUE,
        ifelse(
          m_dt$MatchCount == 2, "B:",
          ifelse(m_dt$MatchCount == 0, "N:",
            ifelse(xor(dt$should_reverse, m_dt$strand == "-"),
              "T:", "U:"
            )
          )
        ),
        "Q:"
      )
      new_levels <- vector_to_combination(
        c("T:", "U:", "B:", "N:", "Q:"),
        levels(dt[[component_col]])
      )
      dt[[component_col]] <- paste0(dt$transcript_bias_label, dt[[component_col]])
      dt[[component_col]] <- factor(dt[[component_col]], levels = new_levels)
    }
  }

  dt.summary <- dt[, .N, by = c(samp_col, component_col)]
  mat <- as.data.frame(data.table::dcast(dt.summary,
    formula = as.formula(paste(samp_col, "~", component_col)),
    fill = 0, value.var = "N", drop = FALSE
  ))

  rownames(mat) <- mat[, 1]
  mat <- mat[, -1]
  mat
}

vector_to_combination <- function(..., c_string = "") {
  expand.grid(
    ...,
    stringsAsFactors = FALSE
  ) %>%
    apply(1, paste0, collapse = c_string) %>%
    unique()
}

count_repeat <- function(x, sequence, is_downstream = TRUE) {
  if (is_downstream) {
    if (startsWith(sequence, x)) {
      pattern <- paste0(
        "^((",
        x,
        ")+).*$"
      )
      nchar(sub(pattern, "\\1", sequence)) / nchar(x)
    } else {
      return(0L)
    }
  } else {
    if (endsWith(sequence, x)) {
      return(count_repeat(
        Biostrings::reverse(x),
        Biostrings::reverse(sequence)
      ))
    } else {
      return(0L)
    }
  }
}


## Examples:
## ACCAA|TC|TAGCGGC or ACAAC|TC|AAGCGGC
## ACCCA|TATC|TATAGCGGC or ACCCATC|TATC|AAGCGGC
## ACCCA|TAGCCTC|TAGCCTAGCGGC or ACCCAGCCTC|TAGCCTC|AAGCGGC
##
## 1
## count_homology_size("TC", "ACCAA", "TAGCGGC")
## count_homology_size("TC", "ACCAC", "AAGCGGC")
## 3
## count_homology_size("TATC", "ACCCA", "TATAGCGGC")
## count_homology_size("TATC", "ACCCATC", "AAGCGGC")
## 5+
## count_homology_size("TAGCCTC", "ACCCA", "TAGCCTAGCGGC")
## count_homology_size("TAGCCTC", "ACCCAGCCTC", "AAGCGGC")
count_homology_size <- function(x, upstream, downstream) {
  size <- nchar(as.character(x))
  if (size >= 2) {
    x_down <- substr(x, 1, size - 1)
    size_down <- 0
    while (nchar(x_down) > 0) {
      if (startsWith(downstream, x_down)) {
        size_down <- nchar(x_down)
        break()
      } else {
        x_down <- substr(x_down, 1, nchar(x_down) - 1)
      }
    }
    x_up <- substring(x, 2)
    size_up <- 0
    while (nchar(x_up) > 0) {
      if (endsWith(upstream, x_up)) {
        size_up <- nchar(x_up)
        break()
      } else {
        x_up <- substring(x_up, 2)
      }
    }
    max(size_up, size_down)
  } else {
    return(0L)
  }
}

is_all_pyrimidine <- function(x) {
  !isTRUE(nchar(gsub("[CT]", "", x)) > 0)
}

is_all_purine <- function(x) {
  !isTRUE(nchar(gsub("[AG]", "", x)) > 0)
}

purine2pyrimidine <- function(x) {
  x <- chartr("A", "T", x)
  x <- chartr("G", "C", x)
  return(x)
}

utils::globalVariables(
  c(
    "A", "APOBEC_Enriched", "APOBEC_Enrichment",
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
    "tCw_to_T", "tcw", "wGa", "wGa_to_A", "wGa_to_C", "wGa_to_T", "wga",
    "Start", "End", "MutIndex",
    ".SD", "dbs", "dbsMotif", "Reference_Allele", "Tumor_Seq_Allele2",
    "chr", "region", "Start_Position", "Hugo_Symbol",
    "End_Position", "ID_len", "ID_motif", "ID_motif_sp", "ID_type",
    "Variant_Type", "count_homosize", "downstream", "mut_base", "should_reverse",
    "upstream"
  )
)

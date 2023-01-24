#' Save plot
save_my_plot <- function(name,
                         plt,
                         type,
                         h = 12,
                         asp = 1.618,
                         path = plots_dir,
                         format = ".pdf") {
  cowplot::save_plot(
    filename = here::here(
      path,
      stringr::str_glue(type,
                        as.character(name),
                        format,
                        .sep = "_"
      )
    ),
    plot = plt,
    base_height = h,
    base_asp = asp,
    limitsize = FALSE
  )
}


#' Select mrtree resolution
SelectResolution <- function(df) {

  #' Get distance between two resolutions with top ARI score
  GetTopResDiff <- function(dat) {
    tmp.ari <-
      dat |>
      top_n(n = 2, wt = ari) |>
      purrr::pluck(2)
    tmp.res <-
      dat |>
      top_n(n = 2, wt = ari) |>
      purrr::pluck(1)
    tmp.ari <- tmp.ari[1] - tmp.ari[2]
    tmp.res <- tmp.res[1] - tmp.res[2]

    return(c(tmp.ari, tmp.res))
  }

  #' Pick one of two top resolutions with parameters
  PickResParam <- function(dat,
                           ari.dif,       res.dif,
                           ari.thd = .05, res.thd = 0) {
    if (ari.dif < ari.thd & res.dif < res.thd) {
      res <-
        dat |>
        top_n(n = 2, wt = ari) |>
        purrr::pluck(1)
      res <- res[2]
    } else {
      res <-
        dat |>
        top_n(n = 1, wt = ari) |>
        purrr::pluck(1)
    }

    return(res)
  }


  df %<>% as_tibble()
  ein.check <-
    df |>
    top_n(n = 2, wt = ari) |>
    purrr::pluck(2) |>
    purrr::map_lgl(~ .x == 1)

  if (all(ein.check)) {
    df %<>% arrange(-resolution) |> distinct(ari, .keep_all = TRUE)
    c(tmp.ari, tmp.res) %<-% GetTopResDiff(df)
    resK <- PickResParam(df, ari.dif = tmp.ari, res.dif = tmp.res)
  } else {
    df %<>%
      filter(ari != 1) |>
      arrange(-resolution) |>
      distinct(ari, .keep_all = TRUE)
    c(tmp.ari, tmp.res) %<-% GetTopResDiff(df)
    resK <- PickResParam(df, ari.dif = tmp.ari, res.dif = tmp.res)
  }

  return(resK)
}


PCScore <- function(object = srt, PCs = 1:5, score.thresh = 1e-05) {
  pAll <- object[["pca"]]@jackstraw$empirical.p.values
  pAll <- pAll[, PCs, drop = FALSE]
  pAll <- as.data.frame(pAll)
  pAll$Contig <- rownames(x = pAll)
  pAll.l <- reshape2::melt(data = pAll, id.vars = "Contig")
  colnames(x = pAll.l) <- c("Contig", "PC", "Value")
  score.df <- NULL
  for (i in PCs) {
    pc.score <-
      suppressWarnings(
        prop.test(
          x = c(length(
            x = which(x = pAll[,
                               i] <= score.thresh)),
            floor(x = nrow(x = pAll) *
                    score.thresh)), n = c(nrow(pAll), nrow(pAll)))$p.val)
    if (length(x = which(x = pAll[, i] <= score.thresh)) == 0) {
      pc.score <- 1
    }
    if (is.null(x = score.df)) {
      score.df <- data.frame(PC = paste0("PC", i), Score = pc.score)
    } else {
      score.df <- rbind(score.df, data.frame(PC = paste0("PC",
                                                         i), Score = pc.score))
    }
  }
  return(score.df)
}


#' Derive MRTree clustering of Seurat object
DeriveKTree <- function(srt, n.pcs = n_pcs, vseed = reseed, n.cores = n_cores) {
  plan("multiprocess", workers = n.cores)

  srt <- NormalizeData(srt)
  srt <-
    FindVariableFeatures(
      srt,
      selection.method = "vst",
      nfeatures = 3000)

  all.genes <- rownames(srt)
  hvg <- VariableFeatures(srt)
  var_regex <- '^Hla-|^Ig[hjkl]|^Rna|^mt-|^Rp[sl]|^Hb[^(p)]|^Gm'
  hvg <- hvg[str_detect(pattern = var_regex, string = hvg, negate = T)]

  srt <-
    ScaleData(srt,
              features = all.genes,
              vars.to.regress = c("var_regex", "log10GenesPerUMI",
                                  "S.Score", "G2M.Score"))

  srt <-
    RunPCA(srt,
           features = hvg,
           npcs = n.pcs,
           seed.use = vseed,
           verbose = TRUE)

  srt <-
    JackStraw(
      object = srt,
      assay = "RNA",
      reduction = "pca",
      dims = n.pcs,
      num.replicate = 100,
      prop.freq = 0.01,
      maxit = 1000)
  srt <-
    ScoreJackStraw(srt,
                   dims = seq_along(srt[["pca"]]@stdev))
  test_pc <-
    PCScore(object = srt,
            PCs = seq_along(srt[["pca"]]@stdev),
            score.thresh = 1e-05)
  selected_pcs <-
    seq_along(
      srt[["pca"]]@stdev
    )[test_pc$Score <= 1e-03 &
        srt[["pca"]]@stdev > quantile(srt[["pca"]]@stdev, .25)]
  selected_pcs

  srt <-
    srt |>
    FindNeighbors(
      dims = selected_pcs,
      k.param = 15,
      annoy.metric = "euclidean",
      n.trees = 100,
      verbose = FALSE) |>
    RunUMAP(
      dims = selected_pcs,
      reduction.name = "umap",
      reduction.key = "UMAP_",
      return.model = FALSE,
      umap.method = "umap-learn",
      densmap = TRUE,
      dens.lambda = 1L,
      dens.frac = 0.3,
      n.epochs = 1000L,
      n.neighbors = 15L,
      min.dist = 0.01,
      spread = 2L,
      metric = "correlation",
      init = "pca",
      seed.use = vseed,
      verbose = FALSE)

  resolutions <-
    modularity_event_sampling(
      A = srt@graphs$RNA_snn,
      n.res = 30,
      gamma.min = 0.2,
      gamma.max = 2.50001
    ) # sample based on the similarity matrix

  srt <- FindClusters(
    srt, algorithm = 4, method = "igraph",
    resolution = resolutions, random.seed = vseed,
    verbose = FALSE)

  out <-  mrtree(
    srt,
    prefix = 'RNA_snn_res.',
    n.cores = n.cores,
    consensus = FALSE,
    augment.path = FALSE
  )

  # Adjusted Multiresolution Rand Index (AMRI)
  ks.flat <-  apply(
    out$labelmat.flat,
    2,
    FUN = function(x)
      length(unique(x))
  )
  ks.mrtree <-  apply(
    out$labelmat.mrtree,
    2,
    FUN = function(x)
      length(unique(x))
  )
  amri.flat <- sapply(1:ncol(out$labelmat.flat), function(i)
    AMRI(out$labelmat.flat[, i], srt$seurat_clusters)$amri)
  amri.flat <- aggregate(amri.flat, by = list(k = ks.flat), FUN = mean)
  amri.recon <- sapply(1:ncol(out$labelmat.mrtree), function(i)
    AMRI(out$labelmat.mrtree[, i], srt$seurat_clusters)$amri)

  df <- rbind(
    data.frame(
      k = amri.flat$k,
      amri = amri.flat$x,
      method = 'Seurat flat'
    ),
    data.frame(k = ks.mrtree, amri = amri.recon, method = 'MRtree')
  )
  stab.out <- stability_plot(out)
  resK <- SelectResolution(stab.out$df)
  srt$k_tree <- out$labelmat.mrtree[, which.min(
    abs(as.integer(
      str_remove(dimnames(
        out$labelmat.mrtree)[[2]], "K"
      )
    ) - resK)
  )]

  Idents(srt) <- "k_tree"
  if (length(unique(srt$k_tree)) > 1) {
    srt.markers.lr <-
      FindAllMarkers(
        srt,
        assay = "RNA",
        verbose = FALSE,
        random.seed = reseed,
        latent.vars = c("var_regex", "log10GenesPerUMI",
                        "S.Score", "G2M.Score"),
        only.pos = TRUE,
        min.pct = 0.05,
        base = 10,
        logfc.threshold = 0.05,
        test.use = "LR")

    if (length(unique(srt.markers.lr$cluster)) > 1) {
      write_csv(
        srt.markers.lr,
        here(tables_dir,
             sprintf('%s_all-mrk_logreg.csv',
                     unique(srt$orig.ident))))}


    srt.markers.mast <-
      FindAllMarkers(
        srt,
        assay = "RNA",
        verbose = FALSE,
        random.seed = reseed,
        latent.vars = c("var_regex", "log10GenesPerUMI",
                        "S.Score", "G2M.Score"),
        only.pos = TRUE,
        min.pct = 0.05,
        base = 10,
        logfc.threshold = 0.05,
        test.use = "MAST")

    if (length(unique(srt.markers.lr$cluster)) > 1) {
      write_csv(
        srt.markers.lr,
        here(tables_dir,
             sprintf('%s_all-mrk_mast.csv',
                     unique(srt$orig.ident))))}
  }

  plan(sequential)
  return(list(srt, srt.markers.lr, srt.markers.mast))
}


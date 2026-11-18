#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
})
options(timeout = 3600, stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args)) args[[1L]] else "E:/Lung/EGA/PDF/gene_list/TCGA/GEO_LUAD_10_validation"
raw_dir <- file.path(base_dir, "raw_data")
out_dir <- file.path(base_dir, "results")
fig_dir <- file.path(base_dir, "figures")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

risk_coef <- c(COL4A4 = -0.0924491994083898, CPS1 = 0.156287432972139,
               IGF2BP3 = 0.0930443595689527, KIAA0319 = 0.0595514108550072,
               NPAS2 = 0.145550168355183, PTPRH = 0.132588888785409)
risk_genes <- names(risk_coef)
gse_ids <- c("GSE3141", "GSE13213", "GSE14814", "GSE29016", "GSE30219",
             "GSE31210", "GSE37745", "GSE50081", "GSE68465", "GSE72094")
platforms <- c(GSE3141 = "GPL570", GSE13213 = "GPL6480", GSE14814 = "GPL96",
               GSE29016 = "GPL6947", GSE30219 = "GPL570", GSE31210 = "GPL570",
               GSE37745 = "GPL570", GSE50081 = "GPL570", GSE68465 = "GPL96",
               GSE72094 = "GPL570") # GPL15048 retains standard Affymetrix probe IDs.

strip_quote <- function(x) gsub('^"|"$', "", x)
num_value <- function(x) {
  z <- as.character(x); out <- rep(NA_real_, length(z))
  hit <- grepl("-?[0-9]+(\\.[0-9]+)?", z)
  out[hit] <- as.numeric(sub(".*?(-?[0-9]+(\\.[0-9]+)?).*", "\\1", z[hit]))
  out
}

platform_bucket <- function(platform) {
  n <- as.integer(sub("GPL", "", platform))
  if (n < 1000) "GPLnnn" else paste0("GPL", floor(n / 1000), "nnn")
}

get_annotation <- function(platform) {
  dest <- file.path(raw_dir, paste0(platform, ".annot.gz"))
  if (!file.exists(dest) || file.info(dest)$size < 1000) {
    url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/platforms/", platform_bucket(platform), "/",
                  platform, "/annot/", platform, ".annot.gz")
    message("Downloading annotation ", platform)
    download.file(url, dest, mode = "wb", method = "libcurl")
  }
  dest
}

read_gene_map <- function(platform) {
  path <- get_annotation(platform)
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  header <- NULL
  repeat {
    line <- readLines(con, 1L, warn = FALSE)
    if (!length(line)) break
    fields <- strip_quote(strsplit(line, "\t", fixed = TRUE)[[1L]])
    if (length(fields) && fields[1L] == "ID") { header <- fields; break }
  }
  if (is.null(header)) stop("No annotation header for ", platform)
  id_col <- match("ID", header)
  symbol_col <- match("Gene Symbol", header)
  if (is.na(symbol_col)) symbol_col <- grep("gene.?symbol|symbol", header, ignore.case = TRUE)[1L]
  probes <- genes <- character()
  repeat {
    lines <- readLines(con, 5000L, warn = FALSE)
    if (!length(lines)) break
    for (line in lines) {
      f <- strip_quote(strsplit(line, "\t", fixed = TRUE)[[1L]])
      if (length(f) < max(id_col, symbol_col)) next
      syms <- trimws(unlist(strsplit(f[symbol_col], " /// |;|,")))
      syms <- intersect(syms, risk_genes)
      if (length(syms)) { probes <- c(probes, rep(f[id_col], length(syms))); genes <- c(genes, syms) }
    }
  }
  unique(data.frame(probe_id = probes, gene = genes, stringsAsFactors = FALSE))
}

read_expression_subset <- function(path, probe_ids) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, 1L, warn = FALSE)
    if (!length(line)) stop("Expression table not found: ", path)
    if (line == "!series_matrix_table_begin") break
  }
  header <- strip_quote(strsplit(readLines(con, 1L), "\t", fixed = TRUE)[[1L]])
  samples <- header[-1L]; keep <- new.env(hash = TRUE, parent = emptyenv())
  for (p in probe_ids) keep[[p]] <- TRUE
  values <- list()
  repeat {
    lines <- readLines(con, 1000L, warn = FALSE)
    if (!length(lines)) break
    for (line in lines) {
      if (line == "!series_matrix_table_end") break
      f <- strip_quote(strsplit(line, "\t", fixed = TRUE)[[1L]])
      if (!exists(f[1L], keep, inherits = FALSE)) next
      val <- suppressWarnings(as.numeric(f[-1L]))
      if (length(val) == length(samples)) values[[f[1L]]] <- val
    }
    if (any(lines == "!series_matrix_table_end")) break
  }
  mat <- do.call(rbind, values); colnames(mat) <- samples; mat
}

collapse_to_genes <- function(probe_expr, gene_map, sample_ids) {
  sample_ids <- intersect(sample_ids, colnames(probe_expr))
  probe_expr <- probe_expr[, sample_ids, drop = FALSE]
  rows <- list(); selected <- list()
  for (g in risk_genes) {
    p <- intersect(gene_map$probe_id[gene_map$gene == g], rownames(probe_expr))
    if (!length(p)) next
    iq <- apply(probe_expr[p, , drop = FALSE], 1L, IQR, na.rm = TRUE)
    chosen <- p[which.max(iq)]
    rows[[g]] <- probe_expr[chosen, ]
    selected[[g]] <- data.frame(gene = g, probe_id = chosen, IQR = max(iq))
  }
  mat <- do.call(rbind, rows); rownames(mat) <- names(rows)
  attr(mat, "selected") <- do.call(rbind, selected)
  mat
}

parse_gse3141 <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  lines <- character()
  repeat {
    z <- readLines(con, 1L, warn = FALSE)
    if (!length(z) || z == "!series_matrix_table_begin") break
    if (grepl("^!Sample_(geo_accession|characteristics_ch1)", z)) lines <- c(lines, z)
  }
  gsm <- strip_quote(strsplit(lines[grepl("geo_accession", lines)], "\t", fixed = TRUE)[[1L]][-1L])
  ch <- strip_quote(strsplit(lines[grepl("characteristics_ch1", lines)], "\t", fixed = TRUE)[[1L]][-1L])
  data.frame(geo_accession = gsm,
             histology = sub(".*Cell type: ([AS]).*", "\\1", ch),
             time = num_value(sub(".*Surv\\(months\\): ([^;]+);.*", "\\1", ch)),
             event = num_value(sub(".*STATUS\\(0=alive, 1=dead\\): ([01]).*", "\\1", ch)))
}

make_clinical <- function(gse, path) {
  if (gse == "GSE3141") {
    x <- parse_gse3141(path); return(x[x$histology == "A", ])
  }
  x <- read.csv(file.path(raw_dir, paste0(gse, "_clinical_extracted.csv")), check.names = FALSE)
  event_dead <- function(z) as.numeric(tolower(trimws(z)) %in% c("dead", "yes", "1"))
  if (gse == "GSE13213") return(data.frame(geo_accession=x$geo_accession,time=num_value(x$survival_days)/30.4375,event=event_dead(x$status)))
  if (gse == "GSE14814") {x<-x[x$histology_type=="ADC",];return(data.frame(geo_accession=x$geo_accession,time=num_value(x$os_time)*12,event=event_dead(x$os_status)))}
  if (gse == "GSE29016") {x<-x[x$histology=="AC",];return(data.frame(geo_accession=x$geo_accession,time=num_value(x$os)*12,event=num_value(x$osbin)))}
  if (gse == "GSE30219") {x<-x[x$histology=="ADC",];return(data.frame(geo_accession=x$geo_accession,time=num_value(x$follow_up_time_months),event=event_dead(x$status)))}
  if (gse == "GSE31210") {x<-x[x$tissue=="primary lung tumor",];return(data.frame(geo_accession=x$geo_accession,time=num_value(x$days_before_death_censor)/30.4375,event=event_dead(x$death)))}
  if (gse == "GSE37745") {x<-x[x$histology=="adeno",];return(data.frame(geo_accession=x$geo_accession,time=num_value(x$days_to_determined_death_status)/30.4375,event=event_dead(x$dead)))}
  if (gse == "GSE50081") {x<-x[x$histology=="adenocarcinoma",];return(data.frame(geo_accession=x$geo_accession,time=num_value(x$survival_time)*12,event=event_dead(x$status)))}
  if (gse == "GSE68465") {x<-x[x$disease_state=="Lung Adenocarcinoma",];return(data.frame(geo_accession=x$geo_accession,time=num_value(x$months_to_last_contact_or_death),event=event_dead(x$vital_status)))}
  if (gse == "GSE72094") return(data.frame(geo_accession=x$geo_accession,time=num_value(x$survival_time_in_days)/30.4375,event=event_dead(x$vital_status)))
  stop("No clinical rule: ", gse)
}

plot_km <- function(dat, st, outfile) {
  fit <- survfit(Surv(time, event) ~ risk_group, data = dat)
  fmt <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  txt <- paste0("Ref: Low risk\nn=",st$n,", events=",st$events,
                "\nLow=",st$low_n,", High=",st$high_n,
                "\nHigh risk: HR=",sprintf("%.2f",st$HR)," (95% CI ",sprintf("%.2f",st$lower95),"-",sprintf("%.2f",st$upper95),")",
                "\nCox P=",fmt(st$cox_p),"\nLog-rank P=",fmt(st$logrank_p))
  p <- ggsurvplot(fit,data=dat,pval=TRUE,conf.int=FALSE,risk.table=TRUE,risk.table.height=0.28,
                  palette=c("#4DBBD5","#E64B35"),xlab="Time (months)",ylab="Overall survival probability",
                  legend.title="Risk group",legend.labs=c("Low risk","High risk"),
                  title=paste0(st$dataset," LUAD OS by fixed 6-gene risk score"))
  p$plot <- p$plot + theme_bw() + theme(aspect.ratio=1,plot.title=element_text(hjust=.5,size=11),
    axis.text=element_text(size=9,color="black"),axis.title=element_text(size=10),legend.text=element_text(size=8)) +
    annotate("text",x=max(dat$time)*.04,y=.08,label=txt,hjust=0,vjust=0,size=2)
  p$table <- p$table + theme_bw() + theme(axis.text=element_text(size=8,color="black"),axis.title=element_text(size=9))
  pdf(outfile,width=5.8,height=7);print(p);dev.off()
}

maps <- list()
for (p in unique(platforms)) {
  maps[[p]] <- read_gene_map(p)
  write.csv(maps[[p]], file.path(out_dir, paste0(p, "_6gene_probe_map.csv")), row.names = FALSE)
}

summary_rows <- list()
for (gse in gse_ids) {
  message("Processing ", gse)
  matrix_file <- file.path(raw_dir, paste0(gse, "_series_matrix.txt.gz"))
  clin <- make_clinical(gse, matrix_file)
  clin <- clin[complete.cases(clin) & clin$time > 0 & clin$event %in% c(0,1), ]
  p <- platforms[[gse]]; probe <- read_expression_subset(matrix_file, maps[[p]]$probe_id)
  expr <- collapse_to_genes(probe, maps[[p]], clin$geo_accession)
  selected <- attr(expr, "selected")
  write.csv(selected, file.path(out_dir, paste0(gse, "_selected_probes.csv")), row.names = FALSE)
  missing <- setdiff(risk_genes, rownames(expr))
  if (length(missing)) stop(gse, " missing model genes: ", paste(missing, collapse=", "))
  common <- intersect(clin$geo_accession, colnames(expr)); clin <- clin[match(common,clin$geo_accession),]
  expr <- expr[,common,drop=FALSE]
  z <- t(scale(t(expr[risk_genes,,drop=FALSE]))); z[!is.finite(z)] <- 0
  score <- as.numeric(crossprod(risk_coef, z))
  dat <- data.frame(dataset=gse,geo_accession=common,time=clin$time,event=clin$event,risk_score=score)
  cut <- median(dat$risk_score); dat$risk_group <- factor(ifelse(dat$risk_score>cut,"High","Low"),levels=c("Low","High"))
  fit <- coxph(Surv(time,event)~risk_group,data=dat); s<-summary(fit)
  lr <- survdiff(Surv(time,event)~risk_group,data=dat); lp<-pchisq(lr$chisq,length(lr$n)-1,lower.tail=FALSE)
  st <- data.frame(dataset=gse,n=nrow(dat),events=sum(dat$event),low_n=sum(dat$risk_group=="Low"),high_n=sum(dat$risk_group=="High"),
    HR=s$conf.int[1,"exp(coef)"],lower95=s$conf.int[1,"lower .95"],upper95=s$conf.int[1,"upper .95"],
    cox_p=s$coefficients[1,"Pr(>|z|)"],logrank_p=lp,median_cutoff=cut,genes_used=paste(risk_genes,collapse=";"))
  summary_rows[[gse]] <- st
  write.csv(dat,file.path(out_dir,paste0(gse,"_fixed_6gene_risk_score_samples.csv")),row.names=FALSE)
  plot_km(dat,st,file.path(fig_dir,paste0("KM_",gse,"_LUAD_OS_fixed_6gene_score.pdf")))
  rm(probe,expr,z);gc(FALSE)
}

tab <- do.call(rbind,summary_rows);rownames(tab)<-NULL;tab$FDR_logrank<-p.adjust(tab$logrank_p,"BH")
write.csv(tab,file.path(out_dir,"GEO_LUAD_10_fixed_6gene_validation_summary.csv"),row.names=FALSE)

forest <- tab; forest$dataset <- factor(forest$dataset,levels=rev(forest$dataset))
pf <- ggplot(forest,aes(HR,dataset,color=logrank_p<.05)) + geom_vline(xintercept=1,linetype="dashed",color="grey40") +
  geom_errorbar(aes(xmin=lower95,xmax=upper95),orientation="y",width=.2) + geom_point(size=2.6) + scale_x_log10() +
  scale_color_manual(values=c(`TRUE`="#E64B35",`FALSE`="#7A7A7A"),guide="none") + theme_classic(base_size=11) +
  labs(title="Fixed TCGA-LUAD 6-gene score in 10 GEO-LUAD cohorts",x="High vs low risk HR (log scale)",y=NULL) +
  theme(plot.title=element_text(hjust=.5,face="bold"),axis.text=element_text(color="black"))
ggsave(file.path(fig_dir,"GEO_LUAD_10_fixed_6gene_HR_forest.pdf"),pf,width=7,height=5.8)
ggsave(file.path(fig_dir,"GEO_LUAD_10_fixed_6gene_HR_forest.png"),pf,width=7,height=5.8,dpi=300)

writeLines(c("Ten GEO LUAD validation cohorts.","Fixed coefficients were trained in TCGA-LUAD and were not refitted in GEO.",
 "Only adenocarcinoma samples were retained in mixed-histology series.",
 "Within each GEO cohort, each gene was z-scored; risk score is the fixed coefficient-weighted sum.",
 "High/Low groups use each GEO cohort's median score; endpoint is overall survival."),file.path(base_dir,"README_GEO_LUAD_10.txt"))
message("Done: ",normalizePath(base_dir,winslash="/"))

# find the most overlap between gene.list.org and gene.list.alt using gene synonym
# output is changed gene names of gene.list.alt similar to gene.list.org

synonym.search <- function(gene.list.org, gene.list.alt) {
  data("ncbi")
  ncbi$Synonyms <- as.character(ncbi$Synonyms)
  syn.symbol <- sapply(1:nrow(ncbi), function(i) unlist(strsplit(ncbi$Synonyms[i], "[|]")))
  all.symbol <- sapply(1:nrow(ncbi), function(i) c(as.character(ncbi$Symbol[i]),syn.symbol[[i]]))
  df <- data.frame(value = unlist(all.symbol),
                   index = rep(seq_along(all.symbol), lapply(all.symbol, length)))
  find.idx <- function(val)df$index[match(val, df$value)]
  gene.list.org <- as.character(gene.list.org)
  gene.list.alt <- as.character(gene.list.alt)
  unmatched.org<- gene.list.org[!gene.list.org %in%gene.list.alt]
  unmatched.alt <- gene.list.alt[!gene.list.alt %in% gene.list.org]
  unmatched.org.id <- find.idx(unmatched.org)
  unmatched.alt.id <- find.idx(unmatched.alt)
  df2=data.frame(id=unmatched.org.id[unmatched.org.id%in% unmatched.alt.id],gene.name=unmatched.org[unmatched.org.id%in% unmatched.alt.id])
  df2<- df2[!is.na(df2$id),]
  ind<- unmatched.alt.id[unmatched.alt.id%in% df2$id]
  unmatched.alt.new <- unmatched.alt
  unmatched.alt.new[unmatched.alt.id%in% df2$id] <- as.character(df2[match(ind,   df2$id),2])
  gene.list.alt.new<- gene.list.alt
  gene.list.alt.new[!gene.list.alt %in% gene.list.org] <- unmatched.alt.new
  gene.list.alt.new
}

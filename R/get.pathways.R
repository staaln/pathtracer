#' @title get.pathways
#' @description Finds homo sapiens pathways in a patway database
#' @param pathwaydatabase Currently only "reactome.db" is implemented
#' @param id One of "symbol" or "entrez"
#' @param min.n.genes Minimal number of genes in a pathway
#' @return A list of pathway ids, athway descriptions, and entrez genes in each pathway
#' @examples get.pathways()

get.pathways = function(pathwaydatabase="reactome.db",id=c("symbol","entrez"), min.n.genes=10) {
  library(reactome.db)
  library(hgug4112a.db)
  # Find pathways (n=2192)
  ptwy.ids = unlist(as.list(reactomePATHNAME2ID))
  ptwy.all = as.list(reactomePATHID2EXTID)
  ptwy.hs = ptwy.all[grep("HSA", names(ptwy.all))]
  ptwy.hs.shortID = names(ptwy.hs)
  ptwy.hs.longID = names(ptwy.ids)[match(ptwy.hs.shortID, ptwy.ids)]

  # Remove pathways with too few genes
  keep = which(sapply(ptwy.hs, length) >= min.n.genes)
  ptwy.hs = ptwy.hs[keep]
  ptwy.hs.shortID = ptwy.hs.shortID[keep]
  ptwy.hs.longID = ptwy.hs.longID[keep]

  # Change identifier if required
  if (id[1] == "symbol") {
    # entrez: Entrez ids, names(entrez): gene symbols
    entrez = unlist(as.list(org.Hs.egSYMBOL2EG))
    gid = lapply(ptwy.hs, function(x) sort(names(entrez)[match(x, entrez)]))
  }
    if (id[1] == "entrez") {
    gid = ptwy.hs
  }

  # Remove any duplicattions in gene lists
  gid = lapply(gid, unique)

  # Return result

  list(id=ptwy.hs.shortID, id2=ptwy.hs.longID, genes=gid)

}

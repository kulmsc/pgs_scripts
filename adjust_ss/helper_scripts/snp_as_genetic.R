snp_asGeneticPos <- function(infos.chr, infos.pos, dir = tempdir(), ncores = 1) {

  #assert_package("R.utils")
  #assert_lengths(infos.chr, infos.pos)

  snp_split(infos.chr, function(ind.chr, pos, dir) {
    chr <- attr(ind.chr, "chr")
    basename <- paste0("chr", chr, ".OMNI.interpolated_genetic_map")
    mapfile <- file.path(dir, basename)
    if (!file.exists(mapfile)) {
      url <- paste0("https://github.com/joepickrell/1000-genomes-genetic-maps/",
                    "raw/master/interpolated_OMNI/", basename, ".gz")
      gzfile <- paste0(mapfile, ".gz")
      utils::download.file(url, destfile = gzfile, quiet = TRUE)
      R.utils::gunzip(gzfile)
    }
    map.chr <- bigreadr::fread2(mapfile, showProgress = FALSE)
    ind <- bigutilsr::knn_parallel(as.matrix(map.chr$V2), as.matrix(pos[ind.chr]),
                                   k = 1, ncores = 1)$nn.idx
    map.chr$V3[ind]
  }, combine = "c", pos = infos.pos, dir = dir, ncores = ncores)
}

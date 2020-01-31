f <- list.files(recursive = TRUE)
f <- f[! f %in% f[grep('\\.R$|gel|config|readMe\\.txt', f)]]
invisible(file.remove(f))

d <- list.dirs()
d <- d[grep('MEGA|alignments|tmp|workFiles',d)]
invisible(unlink(d, recursive = TRUE))
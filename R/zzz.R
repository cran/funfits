require(akima)

.First.lib <- function(lib, pkg) {
  library.dynam("funfits", pkg, lib)
}

if(version$minor < "62")
  library.dynam("funfits")


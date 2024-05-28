packages <- c("magrittr", "osfr", "here")

# check, whether library already installed or not - install and load as needed:
apply(as.matrix(packages), MARGIN = 1, FUN = function(x) {
  
  pkg_avail <- nzchar(system.file(package = x))   # check if library is installed on system
  
  if(pkg_avail){
    require(x, character.only = TRUE)             # load the library, if already installed
    
  }else{
    install.packages(x)                           # install the library, if missing
    require(x, character.only = TRUE)             # load after installation
  }
})




# use osfr-package to download PSACR002-data, on IPD-level, pre-cleaned
osfr::osf_retrieve_file("https://osf.io/rwksu") %>%
  osfr::osf_download(path = here("Data/Downloaded/"))

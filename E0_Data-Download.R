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


### Many Labs 1

# download the zipped-directory of Many Labs 1 data, from the OSF
osfr::osf_retrieve_file("https://osf.io/nqg97") %>% 
  osfr::osf_download(path = here("Data/Downloaded Data"), conflicts = "overwrite")
# the data-file is named very generically "Datasets.zip" - it is renamed here to avoid confusion
file.rename(here("Data/Downloaded Data/Datasets.zip"), here("Data/Downloaded Data/ML1_Datasets.zip"))

# unzip the directory and extract only the specific data-file required for anayses
unzip(here("Data/Downloaded Data/ML1_Datasets.zip"),
      files = "Data/CleanedDataset.sav",
      exdir = here("Data/Original Data/ManyLabs1"),
      junkpaths = TRUE)

### PSACR001

# use osfr-package to download PSACR001-data, on IPD-level, pre-cleaned
osfr::osf_retrieve_file("https://osf.io/jecmr") %>%
  osfr::osf_download(path = here("Data/Original Data/PSACR001"))

### PSACR002

# use osfr-package to download PSACR002-data, on IPD-level, pre-cleaned
osfr::osf_retrieve_file("https://osf.io/rwksu") %>%
  osfr::osf_download(path = here("Data/Original Data/PSACR002"))

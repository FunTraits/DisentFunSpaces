################################################################################
#     _     _ _                    _           
#    | |   (_) |                  (_)          
#    | |    _| |__  _ __ __ _ _ __ _  ___  ___ 
#    | |   | | '_ \| '__/ _` | '__| |/ _ \/ __|
#    | |___| | |_) | | | (_| | |  | |  __/\__ \
#    \_____/_|_.__/|_|  \__,_|_|  |_|\___||___/                    
#
################################################################################

packages<-c('ade4','ape','berryFunctions','betapart','biscale','cowplot','data.table',
            'dggridR','dismo','dplyr','exactextractr','fasterize','FD','funrar','future','future.apply','feather','geiger',
            'geomtextpath','ggplot2','ggpubr','gridExtra','grid','lsmeans','missForest',
            'multcomp','mvMORPH','paleotree','paletteer','pals','paran','phylobase',
            'phytools','picante','plotly','plotrix','progressr','psych','purrr','quanteda',
            'ratematrix','RColorBrewer','readr','readxl','rgbif','rnaturalearth','rredlist',
            'sf','shape','stats','taxize','tidyr','traitdata','TPD','vegan','VennDiagram','viridis','wesanderson')

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
lapply(packages, library, character.only = TRUE)

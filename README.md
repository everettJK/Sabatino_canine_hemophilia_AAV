# Sabatino_canine_hemophilia_AAV

# Create data file for all samples aligned to canFam3
/home/opt/R-3.4.0/bin/Rscript AAVengeR/src/AAVengeR.R AAVengeR/configs/Sabatino.config
cp AAVengeR/outputs/canFam3_final/sites.RData ./data/sites.RData

# Create data file for single chain samples aligned to vector for figure generation.
/home/opt/R-3.4.0/bin/Rscript AAVengeR/src/AAVengeR.R AAVengeR/configs/Sabatino.vectors.single.config

find * -type f -size +100M > .gitignore

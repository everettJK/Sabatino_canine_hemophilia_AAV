# Sabatino_canine_hemophilia_AAV
  

# Create data file for all samples aligned to canFam3  
/home/opt/R-3.4.0/bin/Rscript AAVengeR/src/AAVengeR.R AAVengeR/configs/Sabatino.config  
cp AAVengeR/outputs/canFam3_final/sites.RData ./data/sites.RData
  
  
# Create vector specific data files used to highlight integration within vectors.  
/home/opt/R-3.4.0/bin/Rscript AAVengeR/src/AAVengeR.R AAVengeR/configs/Sabatino.vectors.single.config  
/home/opt/R-3.4.0/bin/Rscript AAVengeR/src/AAVengeR.R AAVengeR/configs/Sabatino.vectors.light.config  
/home/opt/R-3.4.0/bin/Rscript AAVengeR/src/AAVengeR.R AAVengeR/configs/Sabatino.vectors.heavy.config  
  
  
# Create a .gitignore file to ensure that large files are not pushed.  
find * -type f -size +100M > .gitignore  

# Genome Browser in shiny for CONREL;

Webpage of the project at [bcglab.cibio.unitn.it/conrel/](http://bcglab.cibio.unitn.it/conrel/). The Genome Browser is containerized in a Singularity image - download of the image is available at the webserver.
All the data are available inside the singularity image.
The genome browser is implemented using a slightly modified version of the "Interactive visulizations for track-based genomic data in R" implemented by [Marlin-Na/TnT](https://github.com/Marlin-Na/TnT) and available at [ddalfovo/TnT](https://github.com/ddalfovo/TnT)

## Changelog

v2 - Working in progress
 - Improved UI: changed the genome/assembly selection mode
 - Added mouse mm10
 - splitted download files into singularity image and data (to choose between desired organism/assembly). Changed the local singularity installation method

v1.1
 - Added hg38 version and download data
 - Fixed bugs
 - Improved UI
 
v1.0
 - Implementation of the genome browser

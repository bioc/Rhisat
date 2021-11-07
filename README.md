# Rhisat
Bioconductor package: an R wrapper for Hisat and AdapterRemoval

The Rhisat package provides an R interface to the hisat short read aligner by Langmead et al. (2015), and the AdapterRemoval package by Schubert, Lindgreen, and Orlando (2016). The Rhisat package allows users to remove adapter sequences from reads, build bowtie2 indexes (.bt2 or .bt2l), and to create hisat alignment files (.sam or .bam).

## Additional Installation Instructions
The package interfaces with the hisat/ hisat-build/ hisat-inspect wrapper scripts provided in the hisat v0.1.6-beta source code. The hisat wrapper script is a Perl script and the hisat-build/hisat-inspect wrapper scripts are Python scripts. Most versions of MacOS and Linux distributions come with a version of Perl and Python pre-installed. On Windows, both Perl and Python do not come pre-installed and they must be downloaded and installed manually. If either Perl or Python are not installed on your system follow the links below to download and install them.

Python: https://www.python.org/downloads/

Perl: https://www.perl.org/get.html

The package also uses samtools to create bam files if it is present on the system. The reason for this is explained under the Bam File Creation heading. However, samtools is completely optional and the package can be used without it. To download samtools follow the link below.

samtools: http://www.htslib.org/download/

## Hisat source package installation
The Rhisat package uses the hisat v0.1.6-beta source code which was obtained from http://ccb.jhu.edu/software/hisat/downloads/hisat-0.1.6-beta-source.zip. The folders doc, example, scripts, and some non-code files were deleted to reduce the package size.

riboseqr_wrapper
================
`riboSeqR <http://bioconductor.org/packages/3.0/bioc/html/riboSeqR.html>`_
integration for `RiboGalaxy <http://ribogalaxy.ucc.ie/>`_ and `Galaxy <http://galaxyproject.org/>`_.

Included tools
--------------
In the order they are run

1. Prepare riboSeqR Input - Prepare riboSeqR format input files from SAM format alignment files.
   The SAM format files should be obtained by aligning Ribo-Seq and RNA-Seq data to the transcriptome.
   (RNA-Seq data is optional but required for Differential translation analysis).

2. Triplet Periodicity - Plot triplet periodicity for different read lengths.

3. Metagene Analysis

4. Plot Ribosome Profile

   [OR]

   Differential Translation Analysis - Get Ribo and RNA-Seq counts with riboSeqR. Perform differential
   translation analysis with baySeq.

Dependencies
------------
Tested on Ubuntu Linux 14.04 LTS, 64-bit. Dependencies should install automatically on Linux 64-bit.

R ``3.1.2``, riboSeqR ``1.0.5``, baySeq ``2.0.50``, rpy2 ``2.3.10``.

How to test
-----------
1. Upload the following test data files from the test-data folder.

   Prepare riboSeqR input (R data file)

   rsem_chlamy236_deNovo.transcripts.fa

2. A workflow with test data is included in this repository. All tools with the exception of "Prepare riboSeqR input"
   can currently be tested using this workflow. Import this workflow into Galaxy.

3. Run workflow

   In Step 1 of the workflow, select "Prepare riboSeqR input (R data file)" as input.

   In Step 2, select rsem_chlamy236_deNovo.transcripts.fa as input.

   Run workflow.


About the test data files
.........................
The included "Prepare riboSeqR input (R data file)" is saved from an R session using sample data included with the
riboSeqR package. The commands used were ::

   library(riboSeqR)

   datadir <- system.file("extdata", package = "riboSeqR")
   chlamyFasta <- paste(datadir, "/rsem_chlamy236_deNovo.transcripts.fa", sep = "")

   fastaCDS <- findCDS(fastaFile = chlamyFasta, startCodon = c("ATG"), stopCodon = c("TAG", "TAA", "TGA"))

   ribofiles <- paste(datadir, "/chlamy236_plus_deNovo_plusOnly_Index", c(17,3,5,7), sep = "")
   rnafiles <- paste(datadir, "/chlamy236_plus_deNovo_plusOnly_Index", c(10,12,14,16), sep = "")

   riboDat <- readRibodata(ribofiles, rnafiles, replicates = c("WT", "WT", "M", "M"))
   save(riboDat, file="Prepare riboSeqR input (R data file)")

rsem_chlamy236_deNovo.transcripts.fa - sample data from the riboSeqR package.

Bugs/Issues?
------------
Please report here https://github.com/vimalkumarvelayudhan/riboseqr_wrapper/issues

Algorithm 
=========

Overview
--------
.. figure:: /_static/images/Algorithm.png
   :width: 55%
   :align: center
   :alt: Algorithm_overview
|
|
gRNA scoring
------------
The gRNAs score is computed from three weights that aim to (1) maximize specificity, (2) minimize cut-to-insert distance and (3) avoid cutting near splice junctions and in 5’ UTRs.

.. figure:: /_static/images/gRNA.png
   :width: 100%
   :align: center
   :alt: gRNA_scoring
| Notes:
| The gRNA specificity score is calculated in three steps:
| (1) Identify all possible off-target hits of a gRNA in the genome with `BWA <https://bio-bwa.sourceforge.net/>`_. 
| (2) Calculate the off-target `MIT guide specificity score <https://www.nature.com/articles/nbt.2647>`_ for each off-target hit.  
| (3) Take the sum of all MIT scores and use formula 100/(100+sum(mitScores)) to calculate the gRNA specificity score.
| 
|
Recoding strategy
-----------------
| Silent mutations are included in the DNA donor to:
| - Prevent recutting the genome after editing. 
| - Facilitate payload insertion when the cut-to-insert distance is inevitably large.  
|
| There are three recoding intesities: "full", "prevent recut", and "none". 
| In "full", the cut-to-insert region is recoded to facilitate payload insertion. The gRNA or split gRNA (disrupted by the payload, creating a protopsacer-half and a PAM-half) are also recoded.
| "Prevent recut" differs from "full" by the lack of recoding in the cut-to-insert region.

.. figure:: /_static/images/recode.png
   :width: 100%
   :align: center
   :alt: Recode_strategy
      
| Notes
| - The Cutting Frequency Determination (CFD) score was created by `Doench et al. <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4744125/>`_  to calculate the off-target potential of sgRNA:DNA interaction.
|
| Region definitions:
|

.. figure:: /_static/images/region.png
   :width: 100%
   :align: center
   :alt: region_definition
|
|
DNA donor processing strategy
-----------------------------
| After recoding, the DNA donors are further processed, in a type-specific way.
| There are two types of DNA donors:
| - Double-stranded DNA (dsDNA) 
| - Single-stranded oligonucleotides (ssODN)

.. figure:: /_static/images/donor.png
   :width: 100%
   :align: center
   :alt: Donor_strategy

   
.. autosummary::
   :toctree: generated

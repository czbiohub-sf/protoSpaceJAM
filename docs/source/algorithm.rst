Algorithm 
=========

Key concepts
------------
.. figure:: /_static/images/keyConcepts.png
   :width: 100%
   :align: center
   :alt: key_concepts
|
|
Tunable parameters
------------------
.. figure:: /_static/images/tunable_parameters.png
   :width: 100%
   :align: center
   :alt: key_concepts
|
|
gRNA scoring
------------
To rank all candidate gRNAs for a possible design, protoSpaceJAM uses a composite ranking score that weighs (1) the on-target specificity of each candidate, (2) the distance between cut and insertion sites, and (3) the position of the gRNA with respect to important gene expression regulatory sequences, namely 5’ untranslated regions (UTRs) and splice sites  
.. figure:: /_static/images/score.png
   :width: 100%
   :align: center
   :alt: gRNA_scoring

.. figure:: /_static/images/gRNA.png
   :width: 100%
   :align: center
   :alt: gRNA_scoring
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

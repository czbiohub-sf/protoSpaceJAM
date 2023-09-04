Algorithm 
=========
|
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
To rank all candidate gRNAs for a possible design, protoSpaceJAM uses a composite ranking score that weighs (1) the on-target specificity of each candidate, (2) the distance between cut and insertion sites, and (3) the position of the gRNA with respect to important gene expression regulatory sequences, namely 5’ untranslated regions (UTRs) and splice sites.  

.. figure:: /_static/images/score.png
   :width: 100%
   :align: center
   :alt: gRNA_scoring
|
| Weights:
.. figure:: /_static/images/gRNA.png
   :width: 100%
   :align: center
   :alt: gRNA_scoring
|
|
Recoding strategy
-----------------
| protoSpaceJAM supports the optional introduction of silent “recoding” mutations in two separate key regions of the HDR donor:  

* The Cas9/gRNA binding site  
  
The Cas9/gRNA binding site may still be present in the homology arm sequences when payload insertion does not destroy the original protospacer. In such cases, knock-in might be impaired because Cas9 might either cut the donor itself during the delivery of reagents in the cell, or re-cut the knock-in allele after DNA repair. This would respectively decrease donor availability or introduce unwanted genomic modifications, negatively impacting knock-in efficiency overall. A well-established practice is therefore to introduce silent mutations to inactivate the gRNA binding site within the HDR donor. protoSpaceJAM uses the Cutting Frequency Determination (CFD) scoring framework established by Doench and colleagues to predict the impact of individual protospacer and PAM mutations on the Cas9/gRNA cutting potential (14). For each gRNA, protoSpaceJAM identifies the minimum number of mutations that would bring the maximal CFD score in the donor sequence below a user-defined threshold (default: 0.03). When recoding within a protein-coding sequence, only silent mutations are used, leveraging maximal sequence divergence between synonymous codons while excluding rare codons. When recoding within a non-coding region, mutations are introduced in up to one of every three bases. No recoding is allowed in the immediate vicinity of splice junctions, to maintain universally conserved sequence motifs.  

* The cut-to-insert region  
  
When having to perform Cas9/gRNA cuts at a distance from the insertion site, introducing silent mutations in the cut-to-insert region prevents the DNA repair tracks from resolving repair before reaching the payload sequence, thereby increasing the rate of payload insertion. protoSpaceJAM supports recoding within the cut-to-insert region, following the rules outlined above for coding and non-coding sequences and excluding recoding at splice junctions. 

|
| There are three recoding intesities: "full", "prevent recut", and "none". 
| In "full", both the Cas9/gRNA binding site and the cut-to-insert region are recoded.
| In "prevent recut", only the Cas9/gRNA binding site is recoded.

| Recoding strategy summary:
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

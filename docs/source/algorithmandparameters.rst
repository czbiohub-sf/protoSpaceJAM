Algorithm and parameters
========================
|
Key concepts
------------
| Key concepts are illustrated in this example of CRISPR knock-in design: a gRNA (blue) targeting the genomic region of interest, and the HDR donor sequence that will template the genomic insertion of a functional payload.  
| The HDR donor contains silent mutations (yellow stripes) that protects the donor from CRISPR-induced DNA strand breaks. Silent mutations (orange stripes) are added to the cut-to-insert region (between the gRNA cutsite and the payload insertion site) in the HDR donor to safeguard the knock-in. Recoded regions are not considered part of the effective homology arms.  

.. figure:: /_static/images/keyConcepts.png
   :width: 100%
   :align: center
   :alt: key_concepts
|
|
Tunable parameters
------------------
The goal of protoSpaceJAM is to streamline the design of both gRNA and donor sequences using a biologically informed set of rules (summarized in the figure below). These are fully described in subsequent sections of this page.  

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
| How weights are calculated:
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
  
| The Cas9/gRNA binding site may still be present in the homology arm sequences when payload insertion does not destroy the original protospacer. In such cases, knock-in might be impaired because Cas9 could either cut the donor itself during the delivery of reagents in the cell, or re-cut the knock-in allele after DNA repair. This would respectively decrease donor availability or introduce unwanted genomic modifications, negatively impacting knock-in efficiency overall. 
|
| A well-established practice is therefore to introduce silent mutations to inactivate the gRNA binding site within the HDR donor. protoSpaceJAM uses the Cutting Frequency Determination (CFD) scoring framework established by Doench and colleagues to predict the impact of individual protospacer and PAM mutations on the Cas9/gRNA cutting potential. For each gRNA, protoSpaceJAM identifies the fewest mutations that would bring the maximal CFD score in the donor sequence below a user-defined threshold (default: 0.03). When recoding within a protein-coding sequence, only silent mutations are used, leveraging maximal sequence divergence between synonymous codons while excluding rare codons. When recoding in a non-coding region, mutations are introduced in up to one of every three bases. No recoding is allowed in the immediate vicinity of splice junctions, to maintain universally conserved sequence motifs.  

* The cut-to-insert region  
  
When performing Cas9/gRNA cuts at a distance from the insertion site, introducing silent mutations in the cut-to-insert region prevents the DNA repair tracks from resolving repair before reaching the payload sequence, thereby increasing the rate of payload insertion. protoSpaceJAM supports recoding within the cut-to-insert region, following the rules outlined above for coding and non-coding sequences and excluding recoding at splice junctions. 

|
| Recoding strategy summary:
.. figure:: /_static/images/recode.png
   :width: 100%
   :align: center
   :alt: Recode_strategy
      
| Notes:
| - There are three recoding intensities: "full", "prevent recut", and "none". 
|   In "full", both the Cas9/gRNA binding site and the cut-to-insert region are recoded.
|   In "prevent recut", only the Cas9/gRNA binding site is recoded.
| - The Cutting Frequency Determination (CFD) score was created by `Doench et al. <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4744125/>`_  to calculate the off-target potential of sgRNA:DNA interaction.
|
|
DNA donor processing strategy
-----------------------------
| A key goal of protoSpaceJAM is to provide the user with “synthesis-ready” donor sequences to streamline the knock-in experimental process. Therefore, the user can choose between two separate donor design modes - dsDNA and ssODN - that use separate donor processing strategies.  

* Double-stranded DNA (dsDNA)

In dsDNA mode, sequence motifs that might be incompatible with commercial dsDNA synthesis are flagged within the final output table. These flags include homopolymeric runs of 10+ As and Ts or 6+ Gs and Cs and extreme GC content (outside of 25-65% GC content globally, or greater than 52% difference in GC content between any given 50-bp stretches).  

* Single-stranded oligonucleotides (ssODN)

For ssODN synthesis, there is typically no restriction in terms of sequence motifs, but rather in overall length. Therefore, the total length of ssODN donor is capped at a user-defined maximum (default: 200 nt).  
ssODN donors require a choice of polarity for the ssDNA strand to be used. The polarity of the ssODN strand is especially important when using gRNAs with a large cut-to-insert distance. By default, protoSpaceJAM automatically selects the polarity of the ssODN strand to be in the favored orientation. To give the user even finer control over the ssODN strand to be used, four other strand selection modes are also available: Cas9/gRNA target vs non-target strand or transcribed vs non-transcribed strand.  

|
| DNA donor processing strategy summary:
.. figure:: /_static/images/donor.png
   :width: 100%
   :align: center
   :alt: Donor_strategy

   
.. autosummary::
   :toctree: generated

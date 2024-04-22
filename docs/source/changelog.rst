Changelog
=========

|

Algorithm
---------
:Date: April 22, 2024 |enhancement 39eae| Users can now omit transcript IDs when specifying custom genomic coordinates. (recoding will not be aware of coding regions and splice junctions)   

:Date: April 21, 2024 |bug fix bae68| gRNA weight are now raised to the power of the scaling factor.

:Date: April 8, 2024 |enhancement b286e| Added an argument to recode only in coding region. 

:Date: April 8, 2024 |enhancement b286e| The scaling factor of gRNA scoring weights can be zero, the corresponding weight will be ignored if the scaling factor is zero.

:Date: April 4, 2024 |bug fix d8898| Fixed an issue preventing homology arm length be less than 200 bp in dsDNA donor mode.

:Date: March 7, 2024 |enhancement 47b8a| Add scaling factors for weights used in gRNA scoring/ranking.

:Date: March 1, 2024 |enhancement fc5b3| Add support for SpCas9-VQR and enAsCas12a.

:Date: June 26, 2023 |enhancement c5c26| Use Bowtie as the default aligner to check unintended PCR products in primer design.

:Date: June 2, 2023 |enhancement 756ad| Add an option to disable penalizing gRNAs that cut in UTRs or near splice junctions.

:Date: June 2, 2023 |enhancement c3828| |enhancement c2717| Output both trimmed and untrimmed dsDNA donors in the results.

:Date: June 2, 2023 |enhancement c2717| |enhancement c3828| Generate a unique name for each gRNA, donor, and trimmed donor.

:Date: May 15, 2023 |enhancement 9aac65| Results now include the strands of gene, gRNA and donor ( + denotes the forward strand and - denotes the reverse strands). 

:Date: April 11, 2023 |enhancement 6b9db| Users can now specify any genomic coordinates in a transcript as the edit site. 

:Date: Feburary 7, 2023 |bug fix 6e2e5| Fixed a bug that causes mutating the same codons for a second time when scanning and fixing emerging cutsites in payload-homology chimeric regions.  

:Date: Feburary 7, 2023 |bug fix 82c61| Fixed a bug where the last one in an array of recoded nucleoties sometimes were not reflected in the final donor

:Date: Feburary 6, 2023 |enhancement 78fd0| Changed the order of recoding in gRNAs and emerging cutsites (in window scan) to: 3'UTR -> codons -> intron -> 5'UTR

:Date: Feburary 4, 2023 |enhancement 7da6a| Changed the order of recoding in gRNAs to: codons -> 3'UTR -> intron -> 5'UTR

:Date: Feburary 1, 2023 |bug fix e20ed| Fixed two bugs in the dsDNA trimming logic: 1) the right arm was not trimming correctly, 2) report synthesis problems remaining after the trimming step.

:Date: January 25, 2023 |enhancement bf85b| Fine-tuned the off-limit range for recoding near junctions. Now avoiding 3bp/6bp from exon/intron side of the exon/intron junction, and 3bp/2bp from intron/exon side of the intron/exon junction.

:Date: January 24, 2023 |bug fix f87f4| Avoid re-mutating by keeping track of mutated bases. Fixed an interference with mutation caused by marking codon-padding sequence with lower-case.

:Date: January 24, 2023 |enhancement f87f4| Scan and undo isolated mutation of "N" in "NGG".

:Date: January 20, 2023 |bug fix 94ea9| Codon mutation won't happen in some gRNAs in recut mode.

:Date: January 20, 2023 |enhancement c2e92| Changed strand names to match the naming convention.

|

Portal
------
:Date: April 22, 2024 |enhancement 11a19| Users can now omit transcript IDs when specifying custom genomic coordinates. (recoding will not be aware of coding regions and splice junctions)   

:Date: April 21, 2024 |enhancement 95c22| Added a warning message when the primer deisgn can't utilize precomputed primers (longer wait times)

:Date: April 21, 2024 |enhancement 95c22| Added a warning message when the PCR amplicon length is shorter than the DNA donor (PCR amplicon could be generated from the DNA donor).

:Date: April 21, 2024 |enhancement 95c22| Added a warning message when the DNA donor length is shorter than the payload (there is no room for homology arms).

:Date: April 16, 2024 |enhancement b0df0| Tunable primer parameters (e.g. amplicon size and Tm range)

:Date: April 15, 2024 |enhancement 480fd| Improved santitization of user input

:Date: April 14, 2024 |enhancement 2bd98| Annotated DNA donor in GenBank format

:Date: March 07, 2024 |enhancement c5d3b| Added buttons to delete individual submissions from the submission list 

:Date: March 07, 2024 |enhancement c5d3b| Added an option to skip the linker sequence

:Date: March 07, 2024 |enhancement c5d3b| User-tunable parameters will remain last used values when users switch between design steps

:Date: March 01, 2024 |bug fix 25946| Fixed a bug where payloads are empty when defined by selecting tag and linker. 

:Date: March 01, 2024 |enhancement 25946| Add a nuclease selection dropdown for SpCas9 and enAsCas12a.

:Date: Janurary 06, 2024 |enhancement 9a879| Add tooltip instruction for terminus offset (e.g. enable insertions at n bp up and downstream of terminus).

:Date: Janurary 06, 2024 |bug fix 9a879| Fixed a parsing issue preventing the use of submission lists csv file (downloaded from a job) in uploading and populating a new submission (this issue didn't affect the example template).

:Date: July 08, 2023 |bug fix aa7eb| Fixed a bug causing editing payload sequence custom genomic coordinate to fail

:Date: July 08, 2023 |enhancement aa7eb| Make automatic the default strand selection mode for ssODN donors.

:Date: July 06, 2023 |enhancement b95b4| Implement a maximum wait time of 5min for each *ad hoc* GenoPrimer design.

:Date: June 26, 2023 |enhancement 74c6d| Add progress indicator for pJAM.

:Date: June 26, 2023 |enhancement 392f8| Add gene ID and gene symbol for primer output

:Date: June 12, 2023 |enhancement 940a1| Add privacy and cookie policy

:Date: May 19, 2023 |enhancement 04401| Updated to the CZ Biohub SF logo, improved helper text in several places.

:Date: May 19, 2023 |bug fix 04401| CSV upload is updated to work with the new columns in the submission list.

:Date: May 15, 2023 |enhancement c7d70| Genotyping primers are fetched from precomputed results, and if not found, are designed on the fly.

:Date: April 15, 2023 |enhancement 6b9db| Changed the interface to accomodate the input of custom genomic coordinates as edit sites. 

:Date: Feburary 16, 2023 |bug fix 004c6| Entry number are now correct when there are 2+ gRNAs for each design. `Associated change: <https://github.com/czbiohub/protoSpaceJAM-portal/commit/68d37db4642fea22d3738ef5c37da3b9331004c6>`_ ProtospaceJAM will read "Entry" from input, and if fails, uses an auto increment

:Date: Feburary 14, 2023 |enhancement 49990| Added a link in the landing page to a Google form to get an invitation code. Complete the name change to "protoSpaceJAM". And several small changes, e.g. 'Launch' -> 'Jam it'. Fixed typos.

:Date: Feburary 9, 2023 |enhancement dced1| Consolidated donor length parameters into one box, and donor recoding parameters into one box.

:Date: Feburary 1, 2023 |enhancement b6b91| Change the default minumn homology arm length (dsDNA) to 200.

:Date: Feburary 1, 2023 |bug fix b6b91| Made "clear example" and "reset button" buttons work correctly, both will reset to the following defaults: Genome: Human, Genes: None, number of gRNA:1, DNA donor type: ssDNA, HA arm length to consider: 500, target strand: non-target strand, recode intensity:full, prioritize recoding in: PAM, minimum homology arm length: 200, enforce maximum donor length: 200, recut cfd threshold: 0.03.

:Date: January 27, 2023 |enhancement f0ad7| Add a maximum limit of 384 entries per submission list.

:Date: January 26, 2023 |enhancement 0c23a| Default changed to "non-target strand" (including the example).

:Date: January 26, 2023 |enhancement 54621| Default changed to "Prioritize recoding in PAM" (including the example).


.. |enhancement 39eae| image:: https://img.shields.io/badge/39eae-enhancement-green
    https://github.com/czbiohub-sf/protoSpaceJAM/commit/01b9c995ece8109cd9204fb0bdaffe672d039eae
.. |enhancement 11a19| image:: https://img.shields.io/badge/11a19-enhancement-green
    https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/736df18677c6c9b8e84ffa418f7aac8db1011a19
.. |enhancement 95c22| image:: https://img.shields.io/badge/95c22-enhancement-green
    https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/6708a930342a255c8fb64eba0b3356111e195c22
.. |bug fix bae68| image:: https://img.shields.io/badge/bae68-bug%20fix-red
    :target: https://github.com/czbiohub-sf/protoSpaceJAM/commit/774961a0824a59e3bb7294b6ed5df8b28f0bae68
.. |enhancement b0df0| image:: https://img.shields.io/badge/b0df0-enhancement-green
    :target: https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/1f5f7ebda71109305a6b0f3c3e0f44a4d15b0df0
.. |enhancement 480fd| image:: https://img.shields.io/badge/480fd-enhancement-green
    :target: https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/944a1779710d5e3333087ac7d94b534fb78480fd
.. |enhancement 2bd98| image:: https://img.shields.io/badge/2bd98-enhancement-green
    :target: https://github.com/czbiohub-sf/protoSpaceJAM/commit/bda4caee590bee33e1d00de9f067698f6382bd98
.. |enhancement b286e| image:: https://img.shields.io/badge/b286e-enhancement-green
    :target: https://github.com/czbiohub-sf/protoSpaceJAM/commit/491a8936eae7760aeb31c5c0cd6c7ad1a50b286e
.. |bug fix d8898| image:: https://img.shields.io/badge/d8898-bug%20fix-red
    :target: https://github.com/czbiohub-sf/protoSpaceJAM/commit/430b678bf7b9411adee1ab7869fbeff6c37d8898
.. |enhancement 47b8a| image:: https://img.shields.io/badge/47b8a-enhancement-green
    :target: https://github.com/czbiohub-sf/protoSpaceJAM/commit/b2027e1dd0073968008b6e55f6efc64f03647b8a
.. |enhancement c5d3b| image:: https://img.shields.io/badge/c5d3b-enhancement-green
    :target: https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/0fecd264e844d4e6903574b6857635288b2c5d3b
.. |enhancement fc5b3| image:: https://img.shields.io/badge/9a879-enhancement-green
    :target: https://github.com/czbiohub-sf/protoSpaceJAM/commit/0b48770f9767a357b78c9c7c251523dba08fc5b3
.. |bug fix 25946| image:: https://img.shields.io/badge/25946-bug%20fix-red
    :target: https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/65fe28e67fcc93e3f9f3d22e671bbb6e18d25946
.. |enhancement 25946| image:: https://img.shields.io/badge/25946-enhancement-green
    :target: https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/65fe28e67fcc93e3f9f3d22e671bbb6e18d25946
.. |bug fix 9a879| image:: https://img.shields.io/badge/9a879-bug%20fix-red
    :target: https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/9c201a0fa5211f42ad5a94699972d21738e9a879
.. |enhancement 9a879| image:: https://img.shields.io/badge/9a879-enhancement-green
    :target: https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/9c201a0fa5211f42ad5a94699972d21738e9a879
.. |bug fix aa7eb| image:: https://img.shields.io/badge/aa7eb-bug%20fix-red
    :target: https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/4a62c8e95684d8283afd5f038ec2c51acbcaa7eb
.. |enhancement aa7eb| image:: https://img.shields.io/badge/aa7eb-enhancement-green
    :target: https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/4a62c8e95684d8283afd5f038ec2c51acbcaa7eb
.. |enhancement b95b4| image:: https://img.shields.io/badge/b95b4-enhancement-green
    :target: https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/2b6f8b1a004049129037773ff1758acaa60b95b4
.. |enhancement c5c26| image:: https://img.shields.io/badge/c5c26-enhancement-green
    :target: https://github.com/czbiohub-sf/GenoPrimer/commit/f63b44bfa67fd7fbd27d11da1a02c794dfdc5c26
.. |enhancement 74c6d| image:: https://img.shields.io/badge/74c6d-enhancement-green
    :target: https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/ef3101aec0e314123ba2cf8ee7bc1c9571574c6d
.. |enhancement 392f8| image:: https://img.shields.io/badge/392f8-enhancement-green
    :target: https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/327481b312b420fccc2c9c5dc0b5982fbd0392f8
.. |enhancement 940a1| image:: https://img.shields.io/badge/940a1-enhancement-green
    :target: https://github.com/czbiohub-sf/protoSpaceJAM-portal/commit/e405e9c998c23af5bce489d46b76f9ee2c9940a1
.. |enhancement c2717| image:: https://img.shields.io/badge/c2717-enhancement-green
    :target: https://github.com/czbiohub/protoSpaceJAM-portal/commit/d3d055816ea35b9936e7937b91889a139e9c2717
.. |enhancement 756ad| image:: https://img.shields.io/badge/756ad-enhancement-green
    :target: https://github.com/czbiohub/protoSpaceJAM/commit/4bb71f3479236704df299a19ed3da731f97756ad
.. |enhancement c3828| image:: https://img.shields.io/badge/c3828-enhancement-green
    :target: https://github.com/czbiohub/protoSpaceJAM/commit/1a24e1ea0251d4a732d5813240742e6420dc3828
.. |enhancement 04401| image:: https://img.shields.io/badge/04401-enhancement-green
    :target: https://github.com/czbiohub/protoSpaceJAM-portal/commit/d388b8d19d7d1468d4463e0b7061dce1af004401
.. |bug fix 04401| image:: https://img.shields.io/badge/04401-bug%20fix-red
    :target: https://github.com/czbiohub/protoSpaceJAM-portal/commit/d388b8d19d7d1468d4463e0b7061dce1af004401
.. |enhancement 9aac65| image:: https://img.shields.io/badge/9aac65-enhancement-green
    :target: https://github.com/czbiohub/protoSpaceJAM/commit/0566a4d2c79d50190e4df1908d374d4bbb9aac65
.. |enhancement c7d70| image:: https://img.shields.io/badge/c7d70-enhancement-green
    :target: https://github.com/czbiohub/protoSpaceJAM-portal/commit/5631fc0dfb6af3d21a48086c3185ebfdd70c7d70
.. |enhancement ec722| image:: https://img.shields.io/badge/ec722-enhancement-green
    :target: https://github.com/czbiohub/protoSpaceJAM-portal/commit/188f96a2a136678df5a08ee4668a9af3ffaec722
.. |enhancement 6b9db| image:: https://img.shields.io/badge/6b9db-enhancement-green
    :target: https://github.com/czbiohub/protoSpaceJAM/commit/8778e69416078ed2f29499d916724aaac126b9db
.. |bug fix 94ea9| image:: https://img.shields.io/badge/94ea9-bug%20fix-red
    :target: https://github.com/czbiohub/protospaceX/commit/3662c9a9b02e958fd3d6f8a94625470b07b94ea9
.. |bug fix f87f4| image:: https://img.shields.io/badge/f87f4-bug%20fix-red
    :target: https://github.com/czbiohub/protospaceX/commit/98ab6e0dc698effa2441542771d7d82abbdf87f4
.. |enhancement f87f4| image:: https://img.shields.io/badge/f87f4-enhancement-green
    :target: https://github.com/czbiohub/protospaceX/commit/98ab6e0dc698effa2441542771d7d82abbdf87f4
.. |enhancement c2e92| image:: https://img.shields.io/badge/c2e92-enhancement-green
    :target: https://github.com/czbiohub/protospaceX/commit/1b7c70cf2eb6ca6ae8f4783b9337d86a5c7c2e92
.. |enhancement f0ad7| image:: https://img.shields.io/badge/f0ad7-enhancement-green
    :target: https://github.com/czbiohub/protospaceX-portal/commit/687f8faab0839d65da990c9bcbc6487100ff0ad7
.. |enhancement bf85b| image:: https://img.shields.io/badge/bf85b-enhancement-green
    :target: https://github.com/czbiohub/protospaceX/commit/820ed9004c8d33136417ff22733d6812571bf85b
.. |enhancement 0c23a| image:: https://img.shields.io/badge/0c23a-enhancement-green
    :target: https://github.com/czbiohub/protospaceX-portal/commit/823eaff78a281fdfd2627dff329974ccee20c23a
.. |enhancement 54621| image:: https://img.shields.io/badge/54621-enhancement-green
    :target: https://github.com/czbiohub/protospaceX-portal/commit/e80b823bbe1f2a95a9afa6655305402203554621
.. |enhancement b6b91| image:: https://img.shields.io/badge/b6b91-enhancement-green
    :target: https://github.com/czbiohub/protospaceX-portal/commit/1fd046d24253d0fdc8d13d5f1ef9c5f6644b6b91
.. |bug fix b6b91| image:: https://img.shields.io/badge/b6b91-bug%20fix-red
    :target: https://github.com/czbiohub/protospaceX-portal/commit/1fd046d24253d0fdc8d13d5f1ef9c5f6644b6b91
.. |bug fix e20ed| image:: https://img.shields.io/badge/e20ed-bug%20fix-red
    :target: https://github.com/czbiohub/protospaceX/commit/67a4e0df5a33b023e2de834039b4fddd416e20ed
.. |enhancement 7da6a| image:: https://img.shields.io/badge/7da6a-enhancement-green
    :target: https://github.com/czbiohub/protospaceX/commit/1b37873b25f1c0f912f2a3c78445933f1887da6a
.. |enhancement 78fd0| image:: https://img.shields.io/badge/78fd0-enhancement-green
    :target: https://github.com/czbiohub/protospaceX/commit/b70c9762a756355697a7643e0c07af70f4f78fd0
.. |bug fix 6e2e5| image:: https://img.shields.io/badge/6e2e5-bug%20fix-red
    :target: https://github.com/czbiohub/protospaceX/commit/d3b5610d73fd75fa89a9948eb80733bf5286e2e5
.. |bug fix 82c61| image:: https://img.shields.io/badge/82c61-bug%20fix-red
    :target: https://github.com/czbiohub/protospaceX/commit/f94f320dbb9fba33fc6927d39bc2db950ce82c61
.. |enhancement dced1| image:: https://img.shields.io/badge/dced1-enhancement-green
    :target: https://github.com/czbiohub/protospaceX-portal/commit/3818cc5f92e26f170251d950cbadad11c04dced1
.. |enhancement 49990| image:: https://img.shields.io/badge/49990-enhancement-green
    :target: https://github.com/czbiohub/protospaceX-portal/commit/b006e6c3280f0ff09a279e35ec93fb7eb3849990
.. |bug fix 004c6| image:: https://img.shields.io/badge/004c6-bug%20fix-red
    :target: https://github.com/czbiohub/protoSpaceJAM-portal/tree/68d37db4642fea22d3738ef5c37da3b9331004c6

.. autosummary::
   :toctree: generated

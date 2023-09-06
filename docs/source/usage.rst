Usage
=====


Basic usage
-----------

protoSpaceJAM's web app proves a step-by-step navigation system, guiding users through each stage of the process - submission configuration, verification, execution, and visualization of CRISPR knock-in designs  


.. figure:: /_static/images/stepper.png
   :width: 90%
   :align: center
   :alt: Stepper navigation

#. Create a submission list by using one or both of these steps: :guilabel:`Build your job` and/or :guilabel:`Upload csv`  
   (Note that these two steps can repeatedly feed the same submission list)

#. Verify the submission list in the :guilabel:`Submission list` stage

#. Execute protoSpaceJAM in the :guilabel:`Jam it!` stage

#. View and download results in the :guilabel:`Results` stage


1. Create a submission list
---------------------------
| There are two ways to create a submission list.
| Each method can be used repeatedly to feed the same list 

Build your Job (interactive)
    | Select a genome, enter ENST IDs and terminus, adjust on-screen parameters, and then click :guilabel:`Add to the submission list`
    
Upload csv (efficient for large submission list)
    | Click :guilabel:`Download example csv` to obtain a template csv file.
    | Customize the csv file and then click :guilabel:`Upload csv` to upload.
    
   
(Optional) Load an example submission list
------------------------------------------

| There are two options to load an example submission list:

    * **Build your Job**, click :guilabel:`Load example`, and then click :guilabel:`Add to the submission list`.
    
    * **upload csv**, click :guilabel:`Download example csv`, and upload by clicking :guilabel:`Upload csv`.


2. Verify the submission list
-----------------------------
| An example of a submission list is shown below. Click :guilabel:`Confirm` to enable launching protoSpaceJAM

.. figure:: /_static/images/SubmissionList.png
   :width: 100%
   :align: left
   :alt: Submission List 

3. Execute protoSpaceJAM
------------------------
| Click :guilabel:`Jam it` to start processing the submission list
.. figure:: /_static/images/launch.png
   :width: 40%
   :align: center
   :alt: launch
   
|
4. View and download results
----------------------------
| The results page should automatically load after the job is completed.

.. figure:: /_static/images/Results.png
   :width: 100%
   :align: left
   :alt: View/Download results


   
.. autosummary::
   :toctree: generated

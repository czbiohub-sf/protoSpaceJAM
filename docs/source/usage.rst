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
There are two ways to create a submission list:

Build your Job
    | Select a genome, enter ENST IDs and terminus, adjust on-screen parameters and then click :guilabel:`Add to the submission list`
    
Upload csv
    | Click :guilabel:`Download example csv` to obtain a template csv file.
    | Customize the csv file and then click :guilabel:`Upload csv` to upload.
    
    | This method is efficient in uploading a larger submission list.
    
Notes
    | Both methods can add to the same submission list repeatedly
    | Each method can be used repeatedly
    | In **Build your Job**, you can click :guilabel:`Add to the submission list` with a slightly changed configuration 
   
(Optional) Load an example submission list
------------------------------------------

| There are two options to load an example submission list:

    * Option 1: in step **Build your Job**, click :guilabel:`Load example`, and then click :guilabel:`Add to the submission list`.
    
    * Option 2: in step **upload csv**, click :guilabel:`Download example csv`, and upload by clicking :guilabel:`Upload csv`.


2. Verify the submission list
-----------------------------
| An example of a submission list is shown below, click :guilabel:`Confirm` to enable launching protoSpaceJAM

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

.. figure:: /_static/images/Results.png
   :width: 100%
   :align: left
   :alt: View/Download results


   
.. autosummary::
   :toctree: generated

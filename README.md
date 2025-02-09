# MESA_Accreting_2.5Msun

**Author**: Kayla Martin, *PhD Candidate, Macquarie University Sydney (Australia)*.

**Date:** 09/02/2025

**Description:** 

This repository contains the MESA code (version [r24.08.1](https://docs.mesastar.org/en/24.08.1/using_mesa.html)) to evolve a 2.5 solar mass star with initial solar-scaled metallicity through the post-AGB phase, while applying time-varying accretion from a circumbinary disk.

The MESA code scripts provided here are modified versions of those retrieved from [Oomen et al. (2019)](https://www.aanda.org/articles/aa/abs/2019/09/aa35853-19/aa35853-19.html), which were originally prepared using an earlier MESA version (r10398). At the time of the creation of this repository, MESA r24.08.1 was the most recent, stable version of MESA (latest instrumentation paper can be viewed at [this link](https://iopscience.iop.org/article/10.3847/1538-4365/acae8d/meta)). Importantly, more recent versions of MESA have undergone a major restructuring from older version releases (e.g., r10398), and contain many new features (see e.g., [changelog](https://docs.mesastar.org/en/24.08.1/changelog.html)). As such, we have completely updated the MESA scripts from [Oomen et al. (2019)](https://www.aanda.org/articles/aa/abs/2019/09/aa35853-19/aa35853-19.html) to be compatible with newer version releases, specifically version r22.11.1 or later. 

When implementing the MESA code scripts provided here, you will need to make a few modifications to a number of files within the standard MESA directory (download instructions found [here](https://docs.mesastar.org/en/latest/installation.html)). Please refer to the modifications_summary.txt file for a detailed description of these modifications, and a step-by-step guide on the usage of this MESA code.

Full documentation for MESA can be viewed [here](https://docs.mesastar.org/). 

If you have any questions or problems regarding the use of the MESA files contained in this repository, please email me at kayla.martin@hdr.mq.edu.au.

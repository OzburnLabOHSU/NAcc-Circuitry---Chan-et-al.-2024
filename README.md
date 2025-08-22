GitHub assocaited with the below manuscript:

DOI: https://doi.org/10.1007/s00213-025-06875-y

Title: Sex differences in nucleus accumbens circuitry engaged with binge-like ethanol drinking

Authors: Amy E. Chan1,2, Justin Q. Anderson1,2, Kolter B. Grigsby3, Bryan E. Jensen2, Andrey E. Ryabinin1, Angela R. Ozburn1,2

Affiliations: 
1. Oregon Health and Science University, Dept. of Behavioral Neuroscience, Portland Alcohol Research Center, Portland, OR, 97239, USA.
2. Veterans Affairs Portland Health Care System, Research and Development Service, Portland, OR, 97239, USA.
3. Washington State University, Dept. of Pharmaceutical Sciences, Spokane, WA, 99202, USA.


R code necessary to recreate analysis, raw data, figures, and supplemental tables present in the manuscript is provided.

Description of contents:

This folder contains analysis scripts, data associated with the whole-brain c-Fos and NAcc circuitry project in C57BL/6J mice. 

Please cite the source paper if the analysis code was use or adapted for other uses. 

Summary:
61 mice (n=15-16/sex/fluid) underwent stereotactic surgery to deliver AAVrg-hSyn-eGFP bilaterally to the nucleus accumbent core. Following recovery, all mice underwent a 4 day drinking in the dark (DID) task, where mice drank either 20% ethanol or water for 2 hours on days 1-3 and 4 hours on day 4. After the 4th day of drinking, periorbital blood samples were collected to determine blood ethanol concentration (BEC). Immediately after the 4th day of drinking, mice were intracardially perfused with PBS followed by 4% paraformaldehyde (PFA) in PBS. Brains were collected and post-fixed overnight in 4% PFA, then switched to PBS and sent to LifeCanvas Technologies for whole-brain clearing and immunolabeling for c-Fos, GFP, and NeuN. SmartAnalytics was used to register brains to the Allen Brain Atlas and quantify whole-brain c-Fos density and c-Fos+GFP colocalization. 

The raw output of SmartAnalytics contains cell density (cells/mm3) data for 1678 Allen Brain Atlas defined regions per mouse, including layer specific information from the cortex, white matter, ventricles, hindbrain structures and separate values for left and right hemispheres. This data was reduced to 426 areas per mouse (213 per hemisphere), eliminating values for white matter and ventricles, and reducing subregion and layer specific information (ex. infralimbic cortex layer I, II/III, V, VIa, and VIb were collapsed into 1 value for this region). Raw c-Fos density data was found to have a right-skewed distribution, so a square-root transformation was applied to all data to better approximate a normal distribution. 

Analysis performed by:
  Amy Chan - amyechan96@gmail.com
  Justin Anderson - andejust@ohsu.edu

Writing performed by:
  Amy Chan - amyechan96@gmail.com
  Angela Ozburn - ozburn@ohsu.edu
  
Experimental design and conceptualization performed by:
  Amy Chan - amyechan96@gmail.com
  Angela Ozburn - ozburn@ohsu.edu
  
Tissue extraction by:
  Amy Chan - amyechan96@gmail.com
  Kolter Grigsby - kolter.grigsby@wsu.edu

Files included here are organized as follows:

/Source_Code: R script used for differential expression and network analysis.
			  
Corresponding authors:

Amy Chan - amyechan96@gmail.com
Angela Ozburn - ozburn@ohsu.edu


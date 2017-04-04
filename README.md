# epithelial_spheroids
Semi-automated analysis of epithelial spheroid polarity in 3D culture

ImageJ / FIJI macros were developed using the ImageJ version 2.0.0-rc-59/1.51k packaged in the FIJI distribution. MATLAB (The MathWorks, Inc.; Natick, MA, United States) scripts were implemented in version R2015b (Version 8.6- 3 Sep 2015). 
ImageJ macros were written to process Zeiss Vision Image ".zvi" or Tagged Image File ".tif" formats with a resolution of 0.1625 µm per pixel (typical size, 256 x 256 pixels for each plane and color). A different scale of the images may require adaptation of pixel values in spheroids.ijm (see troubleshooting; Additional file 1). Usage of the bio-formats FIJI plugin allows adaption to other image formats. 
For training of the decision tree and usage of trained classifiers in MATLAB analysis with Classification Learner App, the "Statistics and Machine Learning Toolbox" is required. Furthermore, in MATLAB the cell2string function is needed to extract the image names from the files list (see troubleshooting; Additional file 1)

Background: 3D culture models are of increasing importance, allowing analysis of cell properties with characteristics close to in vivo conditions. Epithelial spheroid culture is a simple in vitro model, wherein cells develop highly organized structures with apico-basal polarity and lumen. Spheroid culture allows analysis of epithelial defects as observed for example in ciliopathies, rare (mostly) monogenetic diseases often associated with defects in epithelia of the lung, the kidneys, the liver, the brain and / or the eyes. Here, we report an automated method to classify epithelial spheroids in polarity groups. 
Methods: Seeding of Madin-Darby canine kidney (MDCK) epithelial cells on micro-patterned adhesions chips under controlled conditions leads to growth of epithelial spheroids within 3 days (alternatively, analysis can also be performed on spheroids grown in matrigel). Fluorescence staining for apical (gp135) and basolateral (gp58) markers as well as actin and nuclei is used. To date, classification of spheroid characteristics, e.g. polarity and lumen, is performed blinded by 3 individuals. For a well-founded analysis, at least 100 spheroids in three repeat experiments are required. 
Results: To substitute a bias-prone, time-consuming method of analysis, we established an automated classification of spheroid polarity as core indicator of correct epithelial morphogenesis. A combination of fluorescence signal analyses based on ImageJ / FIJI and MATLAB implementations provides a quantitative readout. Signals in equatorial cross-sections of spheroids are transformed, corrected and normalized. Plot of cumulative marker signals versus spheroid radius provides a quantitative, interpretable representation. The classification algorithm evaluates order and distance of the marker-curves in addition to descriptive shape parameters and generates a comprehensive results file. To allow verification of potentially ambiguous group assignments, inconsistencies in characteristics of individual spheroids are reported. 
Conclusions: In the context of epithelial diseases, spheroids constitute an essential in vitro model. Establishment of a computerized, semi-automated classification of spheroid polarity promises a low effort, time saving analysis, allowing unbiased (high-throughput) characterization of epithelial morphology.

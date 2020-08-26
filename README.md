# Rubber-profile-inspection

Demo code of the research paper:

Arso M Vukicevic, Marko Djapan, Petar Todorovic, Milan Erić, Miladin Stefanovic, Ivan Macuzic. Decision Support System for Dimensional Inspection of Extruded Rubber Profiles. IEEE Access. Vol. 7, pp. 112605-112616, DOI: 10.1109/ACCESS.2019.2934561.

Dimensional inspection of extruded rubber profiles represents currently requires a manual measurement and comparison of profiles’ cross section with the corresponding technical drawings. The considered problem  (in automotive industry) is challenging because the rubber is very flexible material, which makes the manual inspection difficult and time consuming. Although some studies aimed to automate this task, their wider application in industry practice remains limited due to the rubber flexibility and variability of shapes that companies have to produce. Starting from the list of requirements acquired from the industry practice, this study proposed the novel semi-automatic solution based on the usage of affordable equipment and dedicated computer vision algorithms. 

The proposed procedure for dimensional inspection of extruded rubber profiles includes the following steps: 
1) image acquisition, 
2) system calibration, 
3) profile segmentation, 
4) landmark registration, and 
5) augmentation of the referent technical drawing over the acquired image. 

The demo was developed by using a single camera, allowing an operator to make the final decision (as requested by industry partner on the project).

![](images/Graphical%20Abstract%20JPG.jpg)

The code was writen in the Matlab 2019b (run the 'DemoCode/arsGomaGUI.m' script).

Instructions and overview of developed features are available in the video 'Demo video/DemoGomaMP4.mp4'.

![](images/Figure%208.jpg)

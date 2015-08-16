##Automatic Text Detection DEMO##

 
 Yinjie Huang
 
 RET Project
 
 University of Central Florida
 
 2014


###CONTENTS:###


- General Information
- Requirements
- Installation
- Usage
- References

==================================

###GENERAL INFORMATION###

This software was written as the demo of RET project Summer 2014. It represents an implementation of
an automatic text detection algorithm, including a complete graphical user interface (GUI). All rights belong to the author.



###REQUIREMENTS###

To run this software, you need to have the following components installed:

- Mathworks MATLAB
- Mathworks Image Processing Toolbox
- VLFeat open source library  http://www.vlfeat.org/



###INSTALLATION###

This software doesn't require any installation. Just drop the files into a folder.



###USAGE###

To run the software, run the file 'Main.m' or type in 'Main' in the MATLAB command window. The script will take care of all the rest and start a graphical user interface. 

The basic usage is as follows:

- Open one image: click "File" and "open" on image.
- First, we need to make sure if the image:
	(1) Text is light and Background is Dark, we assign the value as 0
	(2) Text is Dark and Background is Light, we assign the value as 1
- We can click the buttons "MSER", "Geometric Filtering", "Stroke Width", "Morphing" and "Final Result" to see the result step by step.
- Step "Stroke Width", the "Std/Mean" means the standard variation devided by mean. We set the default as 0.35. If the detection result is not good. Try other number like 0.5 etc.
- Step "Morphing", "Area Size" means the largest connected component area we consider as normal, any larger ones will be rejected and we do not do image morphing on them. The default is 2000. You can try other numbers like 5000.
- You can also get all the results by clicking "Run by One Click" once.



###References###

- Huizhong Chen; Sam S. Tsai; Georg Schroth; David M. Chen; Radek Grzeszczuk; Bernd Girod., "Robust text detection in natural images with edge-enhanced maximally stable extremal regions," Image Processing (ICIP), 2011 18th IEEE International Conference
on, pp.2609--2612, 2011


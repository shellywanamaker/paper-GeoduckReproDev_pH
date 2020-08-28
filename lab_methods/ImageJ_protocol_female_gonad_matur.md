# ImageJ protocol for histology scoring of female gonad maturation


# Set scale
Open one original image and trace the 100um bar with the line tool. Under analyze in the menu bar, click set scale. enter 100 as the known length and um as the units. Then select global and click ok. 

Close image. The scale will be set for all the following images.

# Measurements
For each sample, do the following:

- Open cropped image.

### Total area measurement:  
- CTRL + m to record the total image area measurement

### Non-tissue area measurement: 

- in the menu bar go to image -> adjust -> brightness/contrast 
	- click set and enter the settings specified in "Settings" tab in this google sheet ([20200820_Female_Gonad](https://docs.google.com/spreadsheets/d/1Y_fKujauAGD3cJCmJRYA8NIpa-VI5repjFRGxhO1M1c/edit?usp=sharing)). These brightness/contrast settings will also apply to other measurements for this image.

- in the menu bar go to image -> adjust -> color threshold
	- adjust hue, saturation and brightness to the settings specified in the "Settings" tab in this google sheet ([20200820_Female_Gonad](https://docs.google.com/spreadsheets/d/1Y_fKujauAGD3cJCmJRYA8NIpa-VI5repjFRGxhO1M1c/edit?usp=sharing)). Be sure to use the "no_tissue_hue", "no_tissue_saturation", and "no_tissue_brightness" min and max settings. Note: if the specified setting are not optimal feel free to adjust and record adjustments to optimize. 
		- these settings should be optimal for highlighting the non-tissue area of the image
		- hit CTRL + m to record the non-tissue area measurement

### Egg measurements:
- in the menu bar go to image -> adjust -> color threshold
	- adjust hue, saturation and brightness to the settings specified in the "Settings" tab in this google sheet ([20200820_Female_Gonad](https://docs.google.com/spreadsheets/d/1Y_fKujauAGD3cJCmJRYA8NIpa-VI5repjFRGxhO1M1c/edit?usp=sharing)). Be sure to use the "egg_hue", "egg_saturation", and "egg_brightness" min and max settings. Note: if the specified setting are not optimal feel free to adjust and record adjustments to optimize. 

- in th menu bar to go measure -> Analyze Particles
	- set size range to 1.00 - infinity
	- set circularity to 0.1-1.00
	- change the Show drop down menu to Outlines
	- check off only display results and in situ show
	- go to File -> Save As and select Tiff 
	- go to File -> Open Recent and select the cropped image you are working on
	- go to Image -> Overlay -> Add Image
		- under Image to add select the .tif you just saved
		- set opacity to 50 and click ok 
			- each egg measurement will propagate in the results window
	- manually check egg measurments by visually scanning through each follicle and
		- note the number for eggs that were incorrectly measured in spreadsheet (you will delete these after exporting the measurement data)
		- manually measure any eggs that were incorrectly measured. These can be artifactually large for instance if two eggs were measured as one, or artifactually small for instance if one egg was measured as two.
			- using lasso tool carefully trace the outline of the egg and hit CTRL + m to record measurement

- **save measurement 'Results' as a text file and import into spreadsheet**
	- I have been using a new table for each sample and have one table for cumulative results 	
#### number of eggs: 
- record the number of eggs after removing the eggs you noted as incorrectly measured

#### egg area:
- calculated the sum of all egg measurements

#### mean egg size:
- calculate the average of all egg measurements

#### standard deviation egg size:
- calculate the standard deviation of all egg measurements

### connective tissue area measurement
- calculate the connective tissue area by subtracting the 'no tissue area' and the 'egg area' from the 'total area'



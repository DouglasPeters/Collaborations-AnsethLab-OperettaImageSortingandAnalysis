# AnsethHCSAnalysisPackage

File: HCS_ImageAnalysis_01062021

Functions Called (in order they are called):
--InFocusImage:
	-Uses a gradient-based system to identify the z-plane that is most likely to be in focus. The rationale is that pixel-pixel gradients will be brighter
	 in z-planes with in-focus signal. If the "highest scoring" z-plane is the first or the last z-plane, then the image is thrown out of all subsequent
	 analyses, since the real focal plane is almost always outside of the image stack in these situations.
	-This uses the DAPI channel to find the focal plane.
	-Focal images are created for each channel, and they are max intensity projections of the focal plane and the two adjacent planes. This method tends to
	 more accurately capture focal signal in samples that aren't perfectly flat.
--CellMaskSegmentation:
	-Uses a marker-based watershedding segmentation method to identify the boundaries of Cell Mask signal and split them into "cell objects," which are
	 clusters of one or more cells.
	-Could technically work with any good cytoplasmic/cell membrane marker, but it was optimized for Cell Mask.
--NuclearSegmentation:
	-Uses a marker-based watershedding segmentation method to identify the boundaries of DAPI signal and split touching nuclei into individual objects,
	 when applicable. Also filters out objects that don't lie within "cell objects."
--CellularAnalysis:
	-This is the main feature-extraction function of the script. It extracts the intensity information for each channel, as well as nearest neighbor stats
	 for each nucleus. Importantly, nearest neighbor calculations are limited to neighbors WITHIN THE SAME CELL OBJECT.
--NucTranslocation:
	-This function will only be called if the "TF_Analysis" variable is set to 1.
	-Calculates the ratio of mean intensities between each nucleus and the "nucleus-free" cytoplasm in which the nucleus resides.
--aSMAActivation:
	-This function will only be called if the "aSMA_Analysis" variable is set to 1.
	-Uses a gradient-based system to quantify aSMA "activation" level in each cell object. In theory, cells with more aSMA fibers should yield higher
	 mean gradient values, but this theory hasn't been thoroughly vetted yet.
	**As far as I'm concerned (circa 12.23.2020), the rationale for this module still needs to be tested thoroughly. The activation level appears
	 to correlate really well with mean aSMA intensity, which I'm not sure it always should.**
--WellSort:
	-This function sorts all the information from individual images in the "Results" variable into individual wells. In other words, it combines all
	 fields of view in each well into a single cell of results.
--ClusterSort:
	-This function will only be called if the "SplitClusterResults" variable is set to 1.
	-Uses the "MinClusterSize" variable value to split the results in each well into "in cluster" and "out cluster" subgroups.

Files Saved (some depend on options selected at top of script):
--FILENAME Segmentation.tif
	-This file will only be saved if the "FigSave" variable is set to 1.
	-An image of the figure showing the segmentation of DAPI and CellMask channels, as well as the aSMA gradient analysis.
--AnalysisResults.mat
	-This file contains the "Results" variable. The analysis results are recorded for each image.
--AnalysisResultsSortedByWell.mat
	-This file contains the "Results_Sorted" variable. The analysis results are combined and sorted by well to make subsequent analysis easier.
--AnalysisResultsSortedByWellandOutsideClusters.mat
	-This file will only be saved if the "SplitClusterResults" variable is set to 1.
	-This file contains the "Results_Sorted_OutClusters" variable. The analysis results for each well (from "Results_Sorted") are separated by
	 how many nuclei were detected in that "cell object." The objects in this variable contain FEWER nuclei than the "MinClusterSize" variable value.
--AnalysisResultsSortedByWellandInsideClusters.mat
	-This file will only be saved if the "SplitClusterResults" variable is set to 1.
	-This file contains the "Results_Sorted_InClustesr" variable. The analysis results for each well (from "Results_Sorted") are separated by
	 how many nuclei were detected in that "cell object." The objects in this variable contain EQUAL/MORE nuclei than the "MinClusterSize" variable value.

Variables in each file:
--AnalysisResults.mat / "Results" variable:
	-FileName: The name of the file created by the FileSort script. (r = row, c = column, f = field).
	-Well: The first 6 characters of FileName, indicating the row and column of the image (i.e., the well).
	-FieldofView: The last 2 characters of FileName, indicating which field of view of the image.
	-ZFocus: The in-focus z-plane, as determined by a gradient-based scoring system in the InFocusImage function file.
	-NucNearestNeighborDistance: The distances between each nucleus and the nucleus nearest it (Euclidean by centroid-centroid distance). These distances are 
		WITHIN cell objects ONLY, so if there is only one nucleus detected in a cell object, there will be no nearest neighbor value. Duplicates are also
		removed, so a cell object with two nuclei will only report one nearest neighbor value.
	-TotalNuclei: The total number of nuclei detected in the image.
	-NucleiPerGroup: The number of nuclei detected in each cell object in the image.
	-AdjustedTotalNuclei: The total number of nuclei detected in the image, and adjusted using the median nuclear area value in the image. This theoretically
		helps correct for improperly segmented clumps of nuclei. As far as I know, this value isn't used for any additional analyses in the script unless
		it explicitly says otherwise.
	-AdjustedNucGroupSizes: The number of nuclei detected in each cell object in the image, but using the AdjustedTotalNuclei value.
	-MeanIntensities: The mean intensity values for each cell object in the image. Each column is a different channel (Col #1 = Channel 1, etc.).
	-SumIntensities: The sum (integrated) intensity values for each cell object in the image. Each column is a different channel (Col #1 = Channel 1, etc.).
	-CMNumber: The total number of cell objects detected in the image (based on Cell Mask (CM) signal segmentation).
	-CMAreas: The area (pixels) of each cell object in the image (based on Cell Mask (CM) signal segmentation).
	-CMTotalArea: The sum (integrated) area of all cell objects in the image (based on Cell Mask (CM) signal segmentation).
	-TFNucCytoRatios: The ratio of mean intensities between each nucleus and the cell object in which it resides. Will only be calculated if TF_Analysis = 1.
	-aSMAGradientMeans: The mean gradient value for each cell object in the image. Will only be calculated if aSMA_Analysis = 1.
--AnalysisResultsSortedByWell.mat / "Results_Sorted" variable:
	-Well: Same as in Results variable.
	-NumFields: The number of images that are analyzed for each well (images can be removed from analysis for various QC reasons).
	-TotalNumberofCMObjects: The sum of all "CMNumber" fields (from Results variable) for each well.
	-TotalNumberofNuclei: The sum of all "TotalNuclei" fields (from Results variable) for each well.
	-TotalAdjustedNumberofNuclei: The sum of all "AdjustedTotalNuclei" fields (from Results variable) for each well.
	-CMNumberObjects: The "CMNumber" field (from Results variable) for each image in each well.
	-CMTotalArea: The "CMAreas" field (from Results variable) for each image in each well.
	-NearestNucDistance: The "NucNearestNeighborDistance" field (from Results variable) for each image in each well.
	-TotalNuclei: The "TotalNuclei" field (from Results variable) for each image in each well.
	-NucleiPerGroup: The "NucleiPerGroup" fields (from Results variable) for each image in each well. The field values from each image are appended in order.
	-AdjustedTotalNuclei: The "AdjustedTotalNuclei" field (from Results variable) for each image in each well.
	-AdjustedNucGroupSizes: The "AdjustedNucGroupSizes" fields (from Results variable) for each image in each well. The field values from each image are appended in order.
	-MeanCh1: The 1st column of the "MeanIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-MeanCh2: The 2nd column of the "MeanIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-MeanCh3: The 3rd column of the "MeanIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-MeanCh4: The 4th column of the "MeanIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-SumCh1: The 1st column of the "SumIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-SumCh2: The 2nd column of the "SumIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-SumCh3: The 3rd column of the "SumIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-SumCh4: The 4th column of the "SumIntensities" field (from Results variable) for each image in each well. The field values from each image are appended in order.
	-TFNucCytoRatios: The "TFNucCytoRatios" field (from Results variable) for each image in each well. The field values from each image are appended in order.
		-This field will only be included if TF_Analysis = 1.
	-aSMAGradientMeans: The "aSMAGradientMeans" field (from Results variable) for each image in each well. The field values from each image are appended in order.
		-This field will only be included i aSMA_Analysis = 1.
--AnalysisResultsSortedByWellandOutsideClusters.mat / "Results_Sorted_OutClusters" variable:
	-Well: Same as in Results variable.
	-CMAreas: The "CMAreas" field (from Results variable) for the cell objects in each image in each well that were determined to be outside of clusters.
	-NumberCMObjects: The number of cell objects in each image in each well that were determined to be outside of clusters.
	-PercentofAllCMObjects: The percentage of all cell objects in all images in the well that were determined to be outside of clusters.
	-TotalCMArea: The sum of "CMAreas" field (from Results variable) for the cell objects in each image in each well that were determined to be outside of clusters.
	-PercentofAllCMArea: The percentage of all cell areas in all images in the well that were determined to be outside of clusters.
	-NucleiPerGroup: The "NucleiPerGroup" field (from Results variable) for cell objects in all images in the well that were determined to be outside of clusters.
	-AdjustedNucGroupSizes: The "AdjustedNucGroupSizes" field (from Results variable) for cell objects in all images in each well that were determined to be outside of clusters.
	-TotalNuclei: The sum of "TotalNuclei" field (from Results variable) for the cell objects in all images in the well that were determined to be outside of clusters.
	-PercentofAllNuclei: The percentage of all nuclei in all images in the well that were determined to be outside of clusters.
	-NearestNucDistanceFiltered: The "NucNearestNeighborDistance" field (from Results variable) for cell objects that were determined to be outside of clusters.
	-MeanCh1: The 1st column of the "MeanIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-MeanCh2: The 2nd column of the "MeanIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-MeanCh3: The 3rd column of the "MeanIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-MeanCh4: The 4th column of the "MeanIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-SumCh1: The 1st column of the "SumIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-SumCh2: The 2nd column of the "SumIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-SumCh3: The 3rd column of the "SumIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-SumCh4: The 4th column of the "SumIntensities" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
	-TFNucCytoRatios: The "TFNucCytoRatios" field (from Results variable) for nuclei with cell objects that were determined to be outside of clusters.
		-This field will only be included if TF_Analysis = 1.
	-aSMAGradientMeans: The "aSMAGradientMeans" field (from Results variable) for cell objects in each image and well that were determined to be outside of clusters.
		-This field will only be included if aSMA_Analysis = 1.
--AnalysisResultsSortedByWellandInsideClusters.mat / "Results_Sorted_InClusters" variable:
	-All fields are the same as in the "Results_Sorted_OutClusters" variable, but only cell objects containing a number of nuclei equal or greater than the MinClusterSize
		variable value are included for analysis.

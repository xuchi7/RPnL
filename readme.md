# Readme

It is the official implementation of RPnL.  
The code was previously released on https://xuchi.weebly.com/rpnl.html  
As weebly cannot be stably accessed from China, I migrate the repository to GitHub.  
The link of this repository is https://github.com/xuchi7/RPnL.git  
If you find this code useful, please cite:

[1] C. Xu; L. Zhang; L. Cheng; R. Koch, "Pose Estimation from Line Correspondences: A Complete Analysis and A Series of Solutions," in IEEE Transactions on Pattern Analysis and Machine Intelligence , vol.PP, no.99, pp.1-1
doi: 10.1109/TPAMI.2016.2582162

The paper can be downloaded free (open access) from IEEE Xplore
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7494617&isnumber=435928

```bibtex
@ARTICLE{rpnl2016, 
	author={C. Xu and L. Zhang and L. Cheng and R. Koch}, 
	journal={IEEE Transactions on Pattern Analysis and Machine Intelligence}, 
	title={Pose Estimation from Line Correspondences: A Complete Analysis and A Series of Solutions}, 
	year={2016}, 
	keywords={Camera Pose Estimation;Configuration Analysis;Perspective-3-Line;Perspective-n-Line}, 
	doi={10.1109/TPAMI.2016.2582162}, 
	ISSN={0162-8828}, 
}
```

## Abstract

In this paper we deal with the camera pose estimation problem from a set of 2D/3D line correspondences, which is also known as PnL (Perspective-n-Line) problem. We carry out our study by comparing PnL with the well-studied PnP (Perspective-n-Point) problem, and our contributions are threefold: (1) We provide a complete 3D configuration analysis for P3L, which includes the well-known P3P problem as well as several existing analyses as special cases. (2) By exploring the similarity between PnL and PnP, we propose a new subset-based PnL approach as well as a series of linear-formulation-based PnL approaches inspired by their PnP counterparts. (3) The proposed linear-formulation-based methods can be easily extended to deal with the line and point features simultaneously.

## keywords 
Camera Pose Estimation;Configuration Analysis;Perspective-3-Line;Perspective-n-Line

## Related Papers
[1] Shiqi Li, Chi Xu*, Ming Xie, "A Robust O(n) Solution to the Perspective-n-Point Problem," IEEE Transactions on  Pattern Analysis and Machine Intelligence, vol. 34, no. 7, pp. 1444-1450,  July 2012, doi:10.1109/TPAMI.2012.41

[2] Lilian Zhang, Chi Xu*, Kok-Meng Lee, and Reinhard Koch. "Robust and efficient pose estimation from line correspondences." In Asian Conference on Computer Vision, pp. 217-230. Springer Berlin Heidelberg, 2012.

[3] Shiqi Li, Chi Xu*. A Stable Direct Solution of Perspective-Three-Point Problem [J]. International Journal of Pattern Recognition and Artificial Intelligence, 2011, 25(5): 643-673

## Usage

Please run "main.m" in this folder.

In "main.m", we tested 6 experiments.  
The defined flags for the experiments are as follows:
- a1; % centered (outlier-free)
- a2; % uncentred (outlier-free)
- b1; % nLine = 4
- b2; % nLine = 5
- c1; % centered (with outliers)
- c2; % uncentred (with outliers)

for each experiment, the code does the following:

(1) generate testing data (3d lines, their 2d projectons, and groundtruth), 
	and store it in "testData" folder

(2)	load the testing data, test the compared methods, and store the results 
	(mean & median rotation & translation errors, correct-rate, execution time) 
	in "resData" folder

(3) load the results in "resData" folder, and plot the figures

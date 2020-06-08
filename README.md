# The *REW-ISA* user's guide #
2020/5/12 15:51:07 

## Function description ##
***REW-ISA*: RNA Expression Weighted Iterative Signature Algorithm**

*REW-ISA* is used to find potential local functional blocks (LFB) based on MeRIP-Seq data, where sites are hyper-methylated or hypo-methylated simultaneously across the specific conditions. *REW-ISA* adopts RNA expression levels of each site as weights to make sites of lower expression level less significant.


## Preparation before using *REW-ISA* ##
### 1. Packages that *REW-ISA* depends on ###
**R package:** biclust, isa2, ggplot2, reshape2 and ggpubr.

Install the above R package using the following code:

    install.packages("biclust")
    install.packages("isa2")
    install.packages("ggplot2")
    install.packages("reshape2")
    install.packages("ggpubr")


### 2. Data preparation ###
#### 2.1 Usually ####
As is known, MeRIP-seq data profiles the m<sup>6</sup>A epi-transcriptome by **input and IP data**.Therefore, before using REW-ISA, it is required to prepare the required IP samples as well as input samples. The two sets of data fed into REW-ISA should be **numerical matrices (or dataframe)**.

#### 2.2 Special cases ####
If the data to be analyzed is not MeRIP-Seq sequencing data, **RNA methylation level** and **RNA expression level** are required.



## Instructions for using *REW-ISA* ##
### 1. Optimizing row and column thresholds based on grid search strategy ###

Given the range of row and column thresholds, *REW-ISA* will calculate under each pair of threshold combinations. We recommend that the algorithm be run repeatedly under each pair of thresholds to ensure the stability of the results. In addition, because the grid search takes a long time, it is recommended to set an appropriate number of repeat runs.

To find the optimized parameters for REW-ISA for the given data, use the following code.

***The following REW-ISA parameter optimization algorithm takes a long time, please run carefully!***


    REWISA_result <- REWISA(FPKM_IP = IP_sample, FPKM_INPUT = input_sample,
                            optimization = TRUE, repeat_num = 10, fixed_side = FALSE, 
                    	    thr_col_interval = seq(0.1, 1.5, 0.05), col_step = 0.05,
                     	    thr_row_interval = seq(1, 3, 0.1), row_step = 0.1)

    # supplement:
    # If you want to optimize the parameters, you must set optimization to TRUE.
    # repeat_num indicates the number of times to run REW-ISA repeatedly under each pair of threshold parameter settings.
    # thr_row_interval and thr_col_interval represent the selection range of row and column thresholds.
    # row_step and col_step indicate that the row and column threshold is within the step size of the selection.
    # For the setting of the line threshold, it is recommended to change it in the range of 1 to 3 in steps of 0.1.
	# For the setting of the column threshold, it is recommended to change it in steps of 0.05 within 0.1 to 1.5.

    # Output of REW-ISA parameter optimization result:
	# ASwC: In each repeated calculation, the Average Similarity within Clusters three-dimensional array calculated for each pair of threshold combinations.
	# SDwC: In each repeated calculation, the Standard Deviation within Clusters three-dimensional array calculated for each pair of threshold combinations.
	# LFB_num: In repeated experiments, a three-dimensional array of LFB numbers generated under each pair of threshold combinations
	# ASwC_mean and SDwC_mean: The average value of each repeated calculation result in each pair of threshold combinations.
	# LFB_num_mode: Under the combination of each pair of thresholds, the mode of the number of LFB is generated.
    # find_TR: Optimized row threshold.
    # find_TC: Optimized col threshold.
    # LFB_number: The optimal number of LFB after optimization.

The optimized parameters are calculated automatically according to the final algorithms of ***LFB\_num\_mode***, ***SDwC\_mean*** and ***ASwC\_mean***.

In order to verify the rationality of the selected parameters, **3D graphs** can be made according to ***LFB\_num\_mode***, ***SDwC\_mean*** and ***ASwC\_mean***.



### 2. Run REW-ISA under the optimized threshold to get the concrete LFBs. ###

If you have optimized parameters, you can run REW-ISA, under specific parameters. Use the following code:

    REWISA_bicluster <- REWISA(FPKM_IP = IP_sample, FPKM_INPUT = input_sample, optimization = FALSE, fixed_side = FALSE, 
    		   				   optimal_thr_row = find_TR, optimal_thr_col = find_TC)
 	# supplement:
    # If you do not need to optimize parameters, you must set optimization to FALSE.



## Contact ##
Please contact the maintainer of *REW-ISA* if you have encountered any problems:

**Shutao Chen:** shutao.chen@cumt.edu.cn

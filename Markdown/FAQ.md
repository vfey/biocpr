#### Q1. What is corplot?   


**A:** corplot is an interactive R application developed using the shiny package for plotting correlation heatmaps from gene expression data.


#### Q2: I have Ensembl IDs instead of HGNC symbols will I still be able to use corplot ?

**A:** Yes, corplot can be used for visualizing correlation heatmaps with Ensembl IDs (in place of HGNC symbols). Incase you wish to change it to HGNC symbols, corplot has an inbuilt option to retrieve HGNC symbols for corresponding Ensembl IDs from biomart, you can select this option upon uploading the data and corplot will retireve the missing HGNC symbols for corresponding Ensembl Ids.


#### Q3. Can you explain the tabs in the tool?


**A:** There are 4 mains tabs in the plot,

```
1. DATA TABLE			-	The data from the uploaded file can be viewed in this file.

2. PLOT DATA			-	The data to be plotted can be viewed in this file, the number of genes to be plotted can be increased/decreased using the slide bar (Number of Genes With Highest Variance) in DATA TABLE tab.

3. CORRELATION HEATMAP          -	The heat map is generated and viewed in this tab.

4. CORRELATION MATRIX           -	The correlation matrix generated can be viewed in this tab.

5. INFO - This tab presents a drop down with two tabs - 'FAQ' and 'ABOUT', from which frequently asked questions and information regarding availability, prerequisite packages and running information can be otained.

``` 


#### Q4: What is the input file format?

**A:** The required input file for corplot is tab separated, where columns represent gene expression values and are labelled with sample names and the rows are labelled according to gene symbols.

```
inputFormat.tsv:

HGNC_symbol		Sample1			Sample2			Sample3

Gene1			7.221534140		3.096289636		0.377739835

Gene2			-1.922776449            0.991096096		3.378650875

Gene3			6.474929692		5.982874941		4.702968285
```

#### Q5: Is there any example data provided?

**A:** Yes, an example data file is provided to try out the tool. You can find the file in `corplot/Data/exampleFIle/` folder.


#### Q6: What do colors in the plot indicate?

**A:** The colors in the plot indicate the level of correlation between two genes. Red means maximum correlation value of +1 whereas blue indicates least correlation value of -1.


#### Q7: Is there a size limit for the expression file?

**A:** The default maximum size limit for files are 5 MB in shiny applications, in corplot the maximum size limit has been modified and set to 40 MB. 

#### Q8: The app crashes when clicked on the 'CORRELATION HEATMAP' tab?

**A:** It is difficult to determine the cause of the crash, in this scenario the most probable cause would be error in uploaded expression file format, please refer the input file format explained in question 4. The other causes for the crash could be due to R package dependencies with different versions, in this case it better to re-install packages and clear the environment (in RStudio) and run again.


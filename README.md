# Visualizing-COVID19-Genomic-Sequences

Independent project done through INFO-B 211: Information Infrastructure II as final project. 

**Fall 2021** <br/>
**Programming Language:** Python <br/>
**Description:**
- Python script utilizes Biopython to conduct comparative analysis between 42 complete genomic COVID-19 sequences.

**Background:** <br/>
Retroviruses, viruses that use RNA as its genetic material, convert their single-stranded (ss) RNA genomes into double-stranded (ds) DNAs. These DNA strands then infect and integrate themselves into the host cell. The existing replication processes within the cell further expresses the retroviral genes into the host organism as the viral DNA is replicated and transcribed into mRNA. (Goldman & Landweber, 2016). The novel coronavirus COVID-19 is a ss retrovirus belonging to the genera Betacoronavirus (SARS-CoV-2, SARS-CoV, and MERS CoV). The similarities of COVID-19 to SARS-CoV with a sequence identity similarity of 79.0 percent termed the virus as SARS-CoV-2 (SARS-CoV-2 Sequencing Data, n.d.). The genomes of coronaviruses in this genera are of the largest genomes among all known RNA viruses with GC contents varying approximately from 32 to 43 percent (Mousavizadeh & Ghasemi, 2021). 

Pairwise sequence alignment is used to identify regions of similarity between all characters in two sequences. The amino acid composition of a genome could give some indicator of the proteins and gene expressivity of a sequence by analyzing which amino acids have the highest frequency (Hormoz, 2013). A T-Test is used to determine if there is significant difference between the means of two groups. The Kruskal-Wallis H Test is a one-way ANOVA that tests the distribution of two samples. <br/>

**Steps of development:** 
1. Extracted genomes from folder of txt files.
2. Created two lists to hold txt file names and contents.
3. Metadata was created for each sequence for 1) Number of Base Pairs (length), 2) Guanine-Cytosine (GC) content, and 3) Molecular Weight. A for loop was used to traverse each genome and calculate the three attributes with the values being added to three separate lists. These lists were put into a dictionary and converted into a data frame to be exported as a csv.
4. A two-way table was created so each sequence was scored against the other to calculate their pairwise sequence alignment score. The function pairwise() takes in a genomic sequence and calculates the pairwise alignment between the sequence passed as a parameter and each consequent sequence. Outside of the pairwise() function, a for loop goes through the genome sequence list, runs pairwise(), and adds the returned list into a new list. This list of lists is converted into a data frame and exported as a csv. 
5. The function amino_acid_composition() takes the translated sequence and returns a dictionary of the amino acid frequencies for that sequence. Outside of the amino_acid_composition() function, a for loop traverses each genome sequence and adds the dictionary of amino acid frequencies to a list. 
6. A horizontal bar graph was developed from the GC contents.
7. A grouped bar graph and heatmap were developed to analyze the pairwise sequence alignment two-way table.
8. A grouped bar graph was developed to analyze the amino acid frequency distribution. 
9. A boxplot was developed from frequencies in Threonine to test if there is an outlier with a higher amino acid frequency. 
10. A T-Test and Kruskal-Wallis H Test were performed on two arrays of Threonine frequencies: with and without the outlier.

**Required Files:** <br/>
Genome Dataset -> Folder containing 42 txt files of complete COVID-19 genome sequences
INFO211_Final_Project.py -> Python script that conducts comparative analysis between COVID-19 sequences

**Output Files:** <br/>
- COVID metadata.csv
- Pairwise Alignment Table.csv
- GC Content Graph.png
- Pairwise Alignment Chart.png
- Pairwise Alignment Heatmap.png
- AA Distribution Chart.png
- AA Boxplot.png

**Conclusion:** <br/>
**GC Content-**<br/> The GC contents for coronaviruses range approximately from 32 to 43 percent. This is concurrent with the results from the GC content metadata in *GC Content Graph.png* as all values lied between 33 and 38 percent. In comparison to the whole sample set, Sequences 28, 29, and 42 vary being approximately 4 percent lower than the GC content mean of 37.75 percent.<br/>

**Pairwise Sequence Alignment-**<br/> *Pairwise Alignment Chart.png* shows Sequences 28, 29, and 42 have the lowest alignment scores compared to other sequences. This is concordant to visualization of the GC content for these sequences as well. *Pairwise Alignment Heatmap.png* supports this with a high concentration of sequences aligned within 99 percent and a stark difference in the alignment scores of Sequences 28, 29, and 42 compared to the whole. Since these sequences displayed a visual difference in comparison to the whole sample set, they were extracted from the larger two-way table for closer analysis. From *Pairwise Alignment Comparison.png* Sequence 42 was the least aligned with all other sequences but had the highest alignment score with Sequence 29 at 93.54 percent showing they are similar, but not a perfect match. Sequence 30 was the least aligned with Sequence 42 at 89.84 percent. <br/>

**Amino Acid Frequency Distribution-**<br/> From *AA Distribution Chart.png*, the amino acids most often expressed were Alanine, Cysteine, Glycine, Threonine and Serine. Looking at Threonine, there appeared to be an outlier with a higher frequency. *AA Boxplot.png* shows that this observed value is an outlier since it lies outside the range of the upper whisker. The p-value of 0.78 from a T-Test didn’t indicate statistical significance that the outlier had an impact on the amino acid distribution between sequences. Additionally, the p-value 0.85 from a Kruskal-Wallis H Test did not indicate the outlier had a statistical significance on the distribution. Using normaltest(), the array without the outlier had a slightly more normal distribution than when the outlier was included. This showed the outlier did have some effect on the overall distribution, but still wasn’t statistically significant as there was only a difference of 0.06 between these two tests. <br/>

**References:** <br/>
Goldman, A. D., & Landweber, L. F. (2016). What Is a Genome? PLOS Genetics, 12(7), e1006181. https://doi.org/10.1371/journal.pgen.1006181 <br/>
Hormoz, S. (2013). Amino acid composition of proteins reduces deleterious impact of mutations. Scientific Reports, 3(1), 2919. https://doi.org/10.1038/srep02919 <br/>
Mousavizadeh, L., & Ghasemi, S. (2021). Genotype and phenotype of COVID-19: Their roles in pathogenesis. Journal of Microbiology, Immunology and Infection, 54(2), 159–163. https://doi.org/10.1016/j.jmii.2020.03.022 <br/>
SARS-CoV-2 Sequencing Data: The Devil Is in the Genomic Detail. (n.d.). ASM.Org. Retrieved December 11, 2021, from https://asm.org/Articles/2020/October/SARS-CoV-2-Sequencing-Data-The-Devil-Is-in-the-Gen



info_text_lookup <- function(tab) {
  if (tab == "Main") {
    return(HTML(
      
"The information button will provide context for the various analyses and visualizations presented here.<br><br>
Check it out on various different tabs.
<div style='text-align: center;'>
<img src='main_tab.png' width='750px'>
</div>
"
      
      
      ))
  } else if (tab == "Genes") {
    return(HTML(
"<p>See where a gene is expressed in terminal cell types, progenitors, or its orthology relationship:</p>
<div style='text-align: center;'>
<img src='genes_tab.png' width='750px'>
</div>
<hr>
<h4>TPM vs. TPM plot:</h4>
<p>The transcripts per million (TPM), which is a summary expression value, is plotted here for the orthologs from <i>C. elegans</i> and <i>C. briggsae</i>. As more genes are selected, a matrix of plots is formed that show expression comparisons between orthologs.</p>
<br>
<h4>TPM vs. percent in plot:</h4>
<p>The transcripts per million (TPM) for each gene individually is plotted against the percent of cells that express that gene. The percent in statistic should be viewed with a bit of caution for comparisons between genes, as how deeply a sequencing library is sequenced can bias the observed % in values.</p>
<br>
<h4>Expression bar plot:</h4>
<p>The expression across the terminal cell types or progenitors is shown for each gene individually as transcripts per million (TPM). The cell types are grouped by their major cell class or lineage origin. Using a linear plotting scale, it is straightforward to see large expression differences.</p>
<br>
<h4>Circular lineage plot:</h4>
<p>Within the progenitor tab, you can see the expression of your genes across the development lineage, starting at the 16-cell stage. Here, gene expression is shown as a color scale using the transcripts per million (TPM). The circle follows the Sulston lineage, with the anterior and posterior divisions shown as left and right daughters, respectively.</p>
<br>
<h4>Orthology:</h4>
<p>Upon doing an ortholog search, the app will see if your gene is in an orthogroup. An orthogroup is a collection of orthologous genes. Here, we determined orthologs using OrthoFinder2. If your gene is in one of the identified orthogroups, it will be displayed here with every gene from <i>C. elegans</i> and <i>C. briggsae</i>. Also, it will show whether any of the genes from the two species were classified as being 1:1 orthologs for the purposes of comparison between species. Summary statistics on these determinations are shown below, including our confidence in the call, whether the genes are syntenic, are each others mutual best match by BLAT score search, whether they were each others mutual best match by Smith-Waterman alignment score, and some details on those alignments.</p>
"
      ))
  } else if (tab == "Cell type") {
    return(HTML(
"<p>Access detailed cell type information and what genes are expressed specifically in that cell type (markers).</p>
<hr>
<h4>Cell type selector:</h4>
<p>At top, you will see a selector for the cell class (tissue) and the cell type if you wish to select a different cell to investigate. Below is a panel that shows your cell type selection, and all of the lineages we have associated with this cell type name. Additionally, if the cell type is a terminal cell type, we display additional naming conventions for the cell type (e.g. AMshL and AMshR for AMsh).</p>
<br>
<h4>Cell type page:</h4>
<p>The cell type page snapshot gives summary details about your cell type including the types of genes specifically expressed there, what these genes are, and how similar the cell type is between species. The top left panel shows every 1:1 ortholog and their expression in Transcripts per Million (TPM) from either species, and colored with whether we identified that gene as having specific expression in your cell type (a cell type marker).</p>
<p>In the top right panel, the WormCat tier 1 annotations are displayed for all of the <i>C. elegans</i> markers as both their total count, and their fold-enrichment.</p>
<p>In the middle panel, the identity of markers that are shared between species (black), are only identified in <i>C. elegans</i> (green), or in <i>C. briggsae</i> (blue) are shown, with their expression in both species.</p>
<p>In the bottom set of panels, there are distributions showing various different summary metrics for the full set of cell types in green for <i>C. elegans</i> and blue for <i>C. briggsae</i>. Your chosen cell type is shown as a vertical line in the two colors for the two species, or as a red vertical line for the pairwise comparisons, with a 95% confidence interval around the point estimate. From these distributions you can access how your cell type compares to other in the dataset, how well it was captured in these experiments, and some biological features about the cell type:</p>
<p><b>Cell count:</b> Total number of cells captured and annotated with that cell type identity.</p>
<p><b>Genes detect:</b> Using bootstrapping when calculating the TPM values, we create a 95% confidence interval around our expression value estimate. If the lower confidence interval intersects zero, we consider this gene not confidently detected. Here we show the number of genes detected for every cell type.</p>
<p><b>Marker count:</b> The number of markers identified in each cell type in each species.</p>
<p><b>Ratio of non-1:1 markers:</b> The ratio of the number of markers that have a 1:1 ortholog to no 1:1 ortholog in the other species.</p>
<p><b>Gini coefficient:</b> A statistic to measure inequality, the gini coefficient shows how balanced (0) or unbalanced (1) a transcriptome is. A value of zero means every gene is expressed equally, which a value of one means just one gene is expressed.</p>
<p><b>Mean number of UMI:</b> The mean number of UMI that were detected in each cell of that cell type.</p>
<p><b>Jensen-Shannon distance:</b> A metric that compares the distance between two probability distributions. Here, we use the two species cell type expression distributions for all 1:1 orthologs to calculate the Jensen-Shannon distance. Please refer to the paper for more details.</p>
<p><b>Pearson correlation:</b> The Pearson correlation between the two species cell type transcriptomes.</p>
<p><b>Cosine distance:</b> The cosine distance is another distance metric that is one minus the cosine angle.</p>
<p><b>Between species Differentially Expressed Genes:</b> The number of genes that were determined to be differentially expressed by a Wilcoxon test. A full list of these genes is found on the GitHub page.</p>
<p>The data presented in this cell page summary is also available with the discrete values in the Cell table tab under Tables.</p>
<br>
<h4>Cell type markers:</h4>
<p>If you’re interested in which genes are specifically expressed in your cell type (cell type markers), the Cell type markers tab will give a visual summary of the top hits.</p>
<p>After selecting your cell, you can adjust the total number of top markers in the two species that are displayed. As discussed in the paper, markers are calculated in each species individually using a Wilcoxon test and filtered by the adjusted p-value, log2FC, and expression level (>80 TPM in that cell type). You can adjust by which metric the top markers are sorted by. We default to a combination of specificity (log2fc) and expression level (TPM). On the x-axis, the measure of specificity is displayed (log2fc), which is the log2 fold-change in expression of your cell type versus the rest of the dataset. So for example, for hyp7_C_lineage, the TPM value of the gene in that cell is compared to the TPM of that gene in the rest of the dataset. The different markers are then grouped by their WormCat tier 1 gene function categorization. You can also look at all of the markers in table form through the Marker table in the Table tab.</p>
"
      ))
  } else if (tab == "Tables") {
    return(HTML(
"<p>Browse summary tables including cell type, gene, and marker data.</p>
<hr>
<h4>Cell table:</h4>
<p>The cell table includes summary information about every cell type captured and annotated in this dataset</p>
<p><b>Cell type:</b> The name of the cell type. Can be filtered.</p>
<p><b>Cell class:</b> The name of the major cell class that the cell type was assigned to. Can be filtered.</p>
<p><b>Jensen-Shannon distance:</b> A metric that compares the distance between two probability distributions. Here, we use the two species cell type expression distributions for all 1:1 orthologs to calculate the Jensen-Shannon distance. Please refer to the paper for more details.  A red color indicates a higher distance while a green color indicates more similarity.</p>
<p><b>Pearson correlation:</b> The Pearson correlation between the two species cell type transcriptomes. A red color indicates a higher distance while a green color indicates more similarity.</p>
<p><b>Genes detect:</b> Using bootstrapping when calculating the TPM values, we create a 95% confidence interval around our expression value estimate. If the lower confidence interval intersects zero, we consider this gene not confidently detected. Here we show the number of genes detected for every cell type.</p>
<p><b>Total markers:</b> The number of markers identified in each cell type in each species.</p>
<p><b>1:1 Markers:</b> The number of markers that have a 1:1 ortholog in the other species.</p>
<p><b>Non-1:1 Markers:</b> The number of markers that have no clear 1:1 ortholog in the other species.</p>
<p><b>Cell count:</b> Total number of cells captured and annotated with that cell type identity.</p>
<p><b>Median UMI:</b> The median number of UMI that were detected in each cell of that cell type.</p>
<br>
<h4>Gene table:</h4>
<p>The gene table includes summary information about every gene that was categorized as a 1:1 ortholog in this dataset. Since it is large, it will take a moment to load.</p>
<p><b>Gene Name Short:</b> The common name for each gene used in this dataset.</p>
<p><b>Filter:</b> Some genes are either poorly detected or poorly expressed in the embryo, and thus can provide dubious summary values for the various metics. They can be filtered out here.</p>
<p><b>Jensen-Shannon Distance:</b> A metric that compares the distance between two probability distributions. Here, we use the two species cell type expression distributions for all 1:1 orthologs to calculate the Jensen-Shannon distance. Please refer to the paper for more details. Here the distances are presented for whether the JSD was calculated on the whole dataset, just the terminal cell types, or just the progenitors. A low value (green) indicates similarity, while a high value (red) indicates a divergence.</p>
<p><b>Tau:</b> The tau is a summary metric where a value of zero means completely even expression across the animal and a value of 1 indicates expression in just one cell type. It can give an indication of the breadth of expression across the embryo.</p>
<p><b>TPM Max Terminal:</b> The maximum expression across any terminal cell type.</p>
<p><b>TPM Max Progenitor:</b> The maximum expression across any progenitor.</p>
<p><b>Gene Name WBGene:</b> The WBGene names for the genes.</p>
<p><b>WormCat:</b> Terms from Amy Walker's group that give a functional classification for every gene in the worm genome.</p>
<br>
<h4>Marker table:</h4>
<p>The marker table includes summary information about every cell type marker identified in either species. You can choose which species' cell type markers to display, and which cell type by selecting above.</p>
<p><b>Gene Name Short:</b> The common name for each gene used in this dataset.</p>
<p><b>Cell Type:</b> The name of the cell type for which markers are being shown.</p>
<p><b>Shared:</b> Whether the marker was identified in both species. This column does not account for whether the gene was 'almost' a marker in one species. Red means not shared and green means shared.</p>
<p><b>Orthology:</b> Whether the gene is classified here as a 1:1 ortholog or not.</p>
<p><b>log2 Fold-Change:</b> A measure of specificity of gene expression. The log2 Fold-change is the log2 ratio in expression of your cell type versus the rest of the dataset</p>
<p><b>Adjusted p-value:</b> When calling markers, a Wilcoxon-rank sum test is run to determine whether a gene is differentially expressed in this cell type versus the rest of the dataset. The values are Bonferroni corrected.</p>
<p><b>TPM Cell Type:</b> The gene expression in this cell type.</p>
<p><b>TPM Max Terminal:</b> The maximum expression across any terminal cell type.</p>
<p><b>TPM Max Progenitor:</b> The maximum expression across any progenitor.</p>
<p><b>Tau:</b> The tau is a summary metric where a value of zero means completely even expression across the animal and a value of 1 indicates expression in just one cell type. It can give an indication of the breadth of expression across the embryo.</p>
<p><b>Gene Name WBGene:</b> The WBGene names for the genes.</p>
<p><b>WormCat:</b> Terms from Amy Walker's group that give a functional classification for every gene in the worm genome.</p>

"
      ))
  }
}


gene_vignette_text = HTML(
"<h4>Goals:</h4>
  <ul>
    <li>Discover the orthology relationship of your gene</li>
    <li>Determine which cell types your gene is expressed in</li>
    <li>Find out whether your gene has a conserved expression pattern between species</li>
  </ul>

<p>First, search for your gene by clearing the content under Search Gene, typing your gene name, and select Orthology under the Gene Sub-Tab. You can search for the common name, the cosmid name, or the WBGene name. Knowing whether your gene was duplicated in <i>C. elegans</i> or <i>C. briggsae</i> will help provide context to any changes in gene expression observed between the species.</p>

<p>On the orthology page, if you searched your gene, you will see the ‘Orthogroup’ which contains your gene, its predicted orthology status based on our analyses, and a maximum-likelihood tree of genes in the orthogroup from the elegans supergroup. Any orthology relationship called here (1:1 versus not 1:1) should be viewed as a hypothesis about its evolutionary history. For the purposes of comparing more evolutionary related orthologs, genes that are duplicated in one species, but show evidence of being at the ancestral location and have higher protein sequence conservation were classified as 1:1 genes.</p>

<div style='text-align: center;'>
<img src='orthology_plot_screenshot.png' width='500px'>
</div>

<p>After determining the orthology status of your gene, let’s now look at what cell types your gene and its orthologs are expressed in. In either the Main tab or the Genes tab, search for your gene in the ortholog search.</p>

<p>You can choose to either look in the terminal cell types or the progenitors for your gene’s expression. Let’s start with the terminal gene expression pattern. On the first tab, if your gene has an ortholog in the reciprocal species, you will see the comparison of your gene’s expression pattern in terminal cell types between <i>C. elegans</i> and <i>C. briggsae</i>. Multiple genes from either species can be plotted here depending on the size of your gene’s orthogroup (see <i>ugt-49</i> as an example below). Below we see intestinal expression of <i>ugt-49</i> in both species.</p>

<div style='text-align: center;'>
<img src='tpm_tpm_plot_screenshot.png' width='500px'>
</div>

<p>On the next tab, you will see the TPM vs. percentage of cells expressing this gene as seen on the WormBase gene pages. Please note that these kinds of plots provide a skewed representation of a gene’s expression patterns as the percentage of cells expressing will vary based on how deeply the single-cell sequencing libraries were sequenced. For our case, the C. briggsae libraries were sequenced less deeply, meaning they will systematically show a lower cell percent expressing.</p>

<p>Next, the expression is shown as a bar plot for easy comparison of expression pattern differences between species and orthologs. You might have to scroll to the right to see the full plots. Below is a clear case of expression difference for the gene <i>ceh-53</i>, which has a remarkable difference in expression in the ASG ciliated neuron.</p>

<div style='text-align: center;'>
<img src='bar_plot_screenshot.png' width='500px'>
</div>

<p>Under the progenitors tab, you will additionally see a Circular lineage plot, which displays the shared invariant lineage, with your gene’s expression overlaid. From this view, it is easiest to see whether a change in expression in a particular development has occurred between species. For example, in the below circle plots for <i>aptf-2</i>, <i>Cbr-aptf-2.1</i> and <i>Cbr-aptf-2.2</i>, we can see an expression difference in the MSxap lineages for <i>Cbr-aptf-2.2</i> relative to <i>aptf-2</i>.</p>

<div style='text-align: center;'>
<img src='circle_plot_screenshot.png' width='500px'>
</div>

<p>For all of these visualizations, the plotting scale (log2 or lineage based) or whether the plots shared a scale (individual or shared) can be altered in the navigation panel to the left. For some visualizations, one style of plotting may be advantageous to see clear differences.</p>

<p>Finally, you can dive deeper into a bunch of summary statistics for all of the genes classified as 1:1 orthologs within the Gene table under the Tables tab.</p>

<div style='text-align: center;'>
<img src='gene_table_screenshot.png' width='500px'>
</div>")

cell_type_vignette_text = HTML("
<h4>Goals:</h4>
  <ul>
  <li>Determine how similar your cell type of interest is between species</li>
  <li>Find genes with specific or a high level of expression in your cell type</li>
  </ul>

<p>First look at the Cell type page snapshot by searching for your cell type using either its lineage identity or the name in our dataset (e.g. hyp4_hyp5_hyp6 or ABplaaaapa):</p>

<div style='text-align: center;'>
<img src='cell_page_screenshot.png' width='500px'>
</div>

<p>At top, you will see a selector for the cell class (tissue) and the cell type if you wish to select a different cell to investigate. Below is a panel that shows your cell type selection, and all of the lineages we have associated with this cell type name. Additionally, if the cell type is a terminal cell type, we display additional naming conventions for the cell type (e.g. AMshL and AMshR for AMsh).</p>

<p>Finally, the cell type page snapshot gives summary details about your cell type including the types of genes specifically expressed there, what these genes are, and how similar the cell type is between species.
The top left panel shows every 1:1 ortholog and their expression in Transcripts per Million (TPM) from either species, and colored with whether we identified that gene as having specific expression in your cell type (a cell type marker). Compared to other cell types, it appears that hyp4_hyp5_hyp6 has several genes with high expression (>30,000 TPM) that are shared between species.</p>

<p>In the top right panel, the WormCat tier 1 annotations are displayed for all of the <i>C. elegans</i> markers as both their total count, and their fold-enrichment. For our below example of hyp4_hyp5_hyp6, we can see genes associated with the terms Extracellular Matrix and Lysosome enriched in this cell type.</p>

<p>In the middle panel, the identity of markers that are shared between species (black), are only identified in <i>C. elegans</i> (green), or in <i>C. briggsae</i> (blue) are shown, with their expression in both species.</p>

<p>In the bottom set of panels, there are distributions showing various different summary metrics for the full set of cell types in green for <i>C. elegans</i> and blue for <i>C. briggsae</i>. Your chosen cell type is shown as a vertical line in the two colors for the two species, or as a red vertical line for the pairwise comparisons, with a 95% confidence interval around the point estimate. From these distributions you can access how your cell type compares to other in the dataset, how well it was captured in these experiments, and some biological features about the cell type. Full information on these metrics is available in the Info button.</p>

<p>The data presented in this cell page summary is also available with the discrete values in the Cell table tab under Tables.</p>

<div style='text-align: center;'>
<img src='cell_table_screenshot.png' width='500px'>
</div>

<p>If you’re interested in which genes are specifically expressed in your cell type (cell type markers), the Cell type markers tab will give a visual summary of the top hits.</p>

<div style='text-align: center;'>
<img src='markers_plot_screenshot.png' width='500px'>
</div>

<p>After selecting your cell, you can adjust the total number of top markers in the two species that are displayed. As discussed in the paper, markers are calculated in each species individually using a Wilcoxon test and filtered by the adjusted p-value, log2FC, and expression level (>80 TPM in that cell type). You can adjust by which metric the top markers are sorted by. We default to a combination of specificity (log2fc) and expression level (log2TPM). On the x-axis, the measure of specificity is displayed (log2fc), which is the log2 fold-change in expression of your cell type versus the rest of the dataset. So for example, for hyp7_C_lineage, the TPM value of the gene in that cell is compared to the TPM of that gene in the rest of the dataset. The different markers are then grouped by their WormCat tier 1 gene function categorization. You can also look at all of the markers in table form through the Marker table in the Table tab.</p>

<div style='text-align: center;'>
<img src='marker_table_screenshot.png' width='500px'>
</div>
<p>Within the marker table, you can additionally see if a marker is ‘shared’ between the two species, meaning whether the marker is in both species. You can additionally sort this table by the different metrics to identify ‘top markers.’</p>")

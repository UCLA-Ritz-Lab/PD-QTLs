# PD-QTLs
Project led by Kimberly Paul to assess whether PD risk SNPs associate with various genomic biomarkers

Copy and pasted message from Kim as a placeholder:

I just re-shared access for you to a box folder called "PD QTLs" on the UCLA Health Box/Ritz account. I've uploaded three markdown report files, that basically go over part of what I have done so far. There's repeated code in them, but also new things, I was still in the deciding best course of action/pipeline phase of this project. So I would read over them starting with "PEG EWAS-GWAS", then "meQTL (Caucasian)", then "meQTL_KP", and just skip through the repeat parts in the latter two. Basically, the second one adds QTL association with PD in and the third CPG association with PD in.

But before this, you probably want to read over some of the papers I mentioned so you better understand what the project even is. At the very bottom of the PEG EWAS-GWAS.html there is a table that has references to previous papers that have linked SNP-CPGs (meQTLs) and then you can also click the link to get metabolomics QTLs/ SNP-metabolite (mQTLs). There's a lot of good papers, but I recommend starting with this one: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1041-x

So basically we have two aims here: 
1) connect genetics to our omics arrays (ie describe QTLs we see linking SNPs to both CPGs (this is the code I have here) and metabolites (newly adding this in)
2) Assess whether these QTLs and associated features (CPGs and metabolites) are associated with PD, or enriched in PD patients, can also do co localization analysis (like the paper I recommended)

Probably best at this step once you've had time to at least read the paper and maybe peruse the code to set-up a meeting? I do have some concerns with the results I was seeing thus far (ie more transQTLs than expected especially versus cisQTLs, maybe I set the distance threshold too low or something, but there could also be something QC going on or even an error somewhere here). So we'll want to nail down the analysis to make sure we trust the results.

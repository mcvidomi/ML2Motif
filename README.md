# ML2Motif

# required installations

1. shogun toolbox<br />
    http://www.shogun-toolbox.org/doc/en/3.0.0/installation.html<br />
    install with python interface
2. R
    https://cran.r-project.org/mirrors.html <br />
    install the additional package bioconductor package by enter in R<br />
        source("http://bioconductor.org/biocLite.R")<br />
        biocLite("seqLogo")<br />


# tutorial

Run "demo.py" script.

1. The human splice data set will be downloaded to the folder data<br />
2. The SVM with a weighted degree string kernel is trained (the time of training depends on how much samples you choose from the training data)
3. gPOIM is computed
4. gPOIM is plotted
5. The important motif is extracted and plotted

If you have any questions, please write a mail at marina.vidovic@tu-berlin.de.

Good luck!

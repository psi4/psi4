PDFS=`find ./output -name "paper.pdf" | sort`

pdftk $PDFS cat output output/proceedings.pdf


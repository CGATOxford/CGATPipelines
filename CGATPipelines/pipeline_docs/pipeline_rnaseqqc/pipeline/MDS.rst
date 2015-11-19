========
MDS plot
========
Descriptive text

MDS plot
========

.. report:: RnaseqqcReport.sampleMDS
   :render: r-ggplot
   :statement: aes(x=MD1, y=MD2) +
	       geom_point() +
	       theme_bw() +
	       xlab('') + ylab ('') +
	       theme(
	       axis.text.x=element_text(size=20, angle=90, hjust=1, vjust=0.5),
	       axis.text.y=element_text(size=20))

   Caption text	       

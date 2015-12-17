Expression Distribution plot
============================

.. report:: RnaseqqcReport.ExpressionDistribution
   :render: r-ggplot
   :statement: aes(x=log2rpkm, group=sample_id, colour=sample_id)+
               geom_density()

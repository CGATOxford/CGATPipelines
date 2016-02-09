.. _ExpressionDistribution:

==================================
Distributions of Expression values
==================================

Insert short description of page contents

Summary::
  * Aims of this analysis
  * What inputs/outputs
  * How the results were generated
  * What you should expect
  * Good example
  * Bad example
  * Links to other examples and the reasons that lie behind them

  This should be a text description of what to expect from the figures on this page.  What
  are the take-home messsages, how should the figures be interpreted, etc

The good

.. report:: GoodExample.Tracker
   :render: myRenderer
   :transform: myTransform
   :options: myAesthetics

   Add a comment about the good example.  What represents good data?

The bad

.. report:: BadExample.Tracker
   :render: myRenderer
   :transform: myTransform
   :options: myAesthetics

   Add a comment about the bad example.  What is specifically bad about this example

More bad examples `<http://myBadData.html >`


Expression Distribution plot
============================

.. report:: RnaseqqcReport.ExpressionDistribution
   :render: r-ggplot
   :statement: aes(x=logTPM, group=sample_id, colour=as.factor(sample_id))+
               geom_density()

Commentary
  This will take the form of some active comments.  This will require the report to
  be published so that it is hosted on the CGAT server/ comments on the DISQUS server.

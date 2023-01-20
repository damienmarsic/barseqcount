
Overview
========

``barseqcount`` is a tool for the analysis of DNA barcode sequencing experiments, typically used for biodistribution studies. It has 2 main commands, ``count`` and ``analyze``.
``barseqcount count`` counts barcodes combinations from NGS read files. ``barseqcount analyze`` analyzes the data and displays the results as a collection of plots.

Source code:
 https://github.com/damienmarsic/barseqcount

Python package:
 https://pypi.org/project/barseqcount/

Bug report / feature requests:
 https://github.com/damienmarsic/barseqcount/issues/new/choose


Installation
============

``barseqcount`` can be installed using pip::

    pip install barseqcount

Note that dependencies might need to be installed individually.


Update
======

To update ``barseqcount`` to the latest version, using pip::

   pip install barseqcount -U


Usage
=====


Functions
=========

main()
******

count(args)
***********

analyze(args)
*************

anaconf(fname,args)
*******************

countconf(fname,args)
*********************

find_bc(l,templ,bcr,cl,ctempl,cbcr)
***********************************

fb(l,templ,i,bcr)
*****************

maxmatch(sample,target,probe)
*****************************

override(func)
**************

version()
*********
Displays version and other information::

    python -m barseqcount -v
      Project: barseqcount
      Version: 0.1.2
      Latest update: 2023-01-20
      Author: Damien Marsic, damien.marsic@aliyun.com
      License: GNU General Public v3 (GPLv3)



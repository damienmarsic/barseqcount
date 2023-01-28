
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
::
    barseqcount.py [-h] [-v] {count,analyze} ...

The optional arguments -h/--help and -v/--version allow to show a help message and version information respectively.

Positional arguments ``count`` and ``analyze``  are the two commands that ``barseqcount`` can run. The former processes read files, collects barcode data and saves a barcode distribution file. The latter analyzes the barcode distribution data and displays the results as a collection of plots. 

``barseqcount count`` has optional arguments -c/--configuration_file and -n/--new. The -c argument (followed by a file name) specifies which configuration file to use, or which to create if it does not exist yet. The -n argument allows to ignore an existing configuration file and to create a new one.

``barseqcount analyze`` has the same -c and -n arguments, but also a third one: -f/file_format, allowing to chose a file format to save the plots individually. If the -f argument i not used, then all plots will be saved in a single multipage pdf file.

barseqcount count
*****************


barseqcount analyze
*******************


Functions
=========

Many of the functions used in ``barseqcount`` are also used in other projects and have been iincluded in the `dmbiolib <https://dmbiolib.readthedocs.io/en/latest/dbl-doc.html>`_ package.

main()
******

The ``main()`` function uses ``argparse`` to read and process the command line arguments. 

count(args)
***********
* args: optional arguments following the ``count`` command

| Creates a new configuration file if none exists or if -n/--new argument is present. Otherwise, processes the read file(s) according to instructions in the configuration file. Saves the barcode distribution in a csv file, and a report in a txt file.

analyze(args)
*************
* args: optional arguments following the ``analyze`` command

| Creates a new configuration file if none exists or if -n/--new argument is present. Otherwise, analyzes the data according to instructions in the configuration file. Creates a series of plots and saves resultst in csv files.

anaconf(fname,args)
*******************
* fname: name of the configuration file to be created
* args: arguments

| Creates a configuration file for the ``barseqcount analyze`` program

countconf(fname,args)
*********************
* fname: name of the configuration file to be created
* args: arguments

| Creates a configuration file for the ``barseqcount count`` program

find_bc(l,templ,bcr,cl,ctempl,cbcr)
***********************************
* l: read
* templ: template
* bcr: dictionary containing information about barcode locations and error correction
* cl: compressed read (using compress function from ``dmbiolib``)
* ctempl: compressed template
* cbcr: dictionary containing information about barcode locations based on compressed template

| Identifies all barcodes in a read and perfoems error correction as appropriate.

| Returns a dictionary of barcode positionsa / barcode sequences, a number indicating whether the read was corrected (>0) or not (0), and a list containing error correction counters.

fb(l,templ,i,bcr)
*****************
* l: read (nucleotide sequence)
* templ: template
* i: barcode index
* bcr: dictionary containing information about barcode locations and error correction

| Determines bacode sequence by mapping read sequence to template, using information about barcode locations and error correction.

| Returns barcode sequence.

maxmatch(sample,target,probe)
*****************************
* sample: nucleotide sequence of primer
* target: nucleotide sequence of template
* probe: initial probe size

| Determines largest part of the primer that matches the template.

| Returns (a,x,b,y) where a is the maximum extent of the primer from its right end that matches the template, b is the maximum extent of the primer from its left end that matches the template, x is the template index of sample[-a:], and y is the template index of sample[:b].

override(func)
**************
Allows argparse to handle the -v / --version argument correctly.

version()
*********
Displays version and other information::

    python -m barseqcount -v
      Project: barseqcount
      Version: 0.1.2
      Latest update: 2023-01-20
      Author: Damien Marsic, damien.marsic@aliyun.com
      License: GNU General Public v3 (GPLv3)



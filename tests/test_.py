#!/usr/bin/env python
import pytest
bsc=__import__('barseqcount')

def test_maxmatch():
    assert bsc.maxmatch('gtcaaagcttag','aaacggtcaaagctgtaggcaacatgtcag',5)==(0,0,9,5)

def test_fb():
    assert bsc.fb('tactgcagcttcgtacgggttacct','tactnnnnnttcgtacgggttacct',4,{4:[5,0,9,14]})=='gcagc'

def test_find_bc():
    assert bsc.find_bc('tactgcagcttcgtacgggttacct','tactnnnnnttcgtacgggttacct',{4:[5,0,9,14]},'tactgcagctcgtacgtact','tactnnnnntcgtacgtact',{4:[5,0,9,14]})==({4: 'gcagc'}, 0, [0, 0, 0])

pytest.main()


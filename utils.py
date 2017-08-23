#!/usr/bin/env python
import os
import sys


def format_seconds(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


def welcome_msg(title, author):
	'Print Welcome message'
	print '\n\n'
	print '\t======================================================='
	print '\t\t\t\t\t\t\t\t'
	print '\t\t :-) >>------->',title,'<-------<< (-:	\t'
	print '\t\t\t\t\t\t\t\t'
	print '\t\t\t\t\t\t\t\t'
	print '\t\tThis programe performs the', title,' \t\t| \n\t|\ton a MD trajectory\t\t\t\t'
	print '\t\t\t\t\t\t\t\t|\n', '\t|\tAuthor(s):', author,'\t\t\t\t|\n','\t|\tResearch Unit in Bioinformatics (RUBi)\t\t|\n', '\t|\tRhodes University, 2017\t\t\t\t'
	print '\t\tDistributed under GNU GPL 3.0\t\t\t'
	print '\t\t\t\t\t\t\t\t'
	print '\t\thttps://github.com/michaelglenister/NMA-TASK\t'
	print '\t\t\t\t\t\t\t\t'
	print '\t======================================================='
	print '\n'
	return;
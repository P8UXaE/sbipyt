#!/bin/bash
'''
This file allows you to create automatically the working environment and load the necessary packages to run the biter programm. 

To activate the environment just type 
$ source python3_9venv/bin/activate
'''

# create virtual environment
python3.9 -m venv python3_9venv

# activate virtual environment
# source myenv/bin/activate

# install all the packages 
python setup.py install

# deactivate virtual environment
# deactivate


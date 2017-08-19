Installation
========================================

Platform compatibility
-------------------------------

NMA-TASK is compatible with most platforms which are able to run Python 2.7 and G++


Install system dependencies
-----------------------------

**Ubuntu 16.04:** ::

	sudo apt-get install virtualenvwrapper python-dev g++


Install Python dependencies
--------------------------------

It is recommended to create a Python virtual environment for installing and managing dependencies::

	virtualenv venv
	source venv/bin/activate
	pip install --upgrade pip
	pip install numpy
	pip install matplotlib
	pip install cython
	pip install mdtraj


Download the project
-------------------------------

NMA-TASK can be cloned from it's GitHub repository ::

	git clone https://github.com/RUBi-ZA/NMA-TASK.git
	cd NMA-TASK

Activate the virtual environment you created in the previous step when using NMA-TASK. with::

	source venv/bin/activate

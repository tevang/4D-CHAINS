Brief Description
=================

4D-CHAINS is a software for fully automated protein backbone N-H and sidechain aliphatic C & H chemical shift assignment from 2 NMR spectra: a 4D TOCSY and
a 4D NOESY. 


License
============

4D-CHAINS software for protein NMR assignment is a property of Thomas Evangelidis and Konstantinos Tripsianes and is free for NON-COMMERCIAL USAGE. The code is licensed under the Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY- NC-ND 4.0). You are free to:

* Share - copy and redistribute the material in any medium or format.
* The licensor cannot revoke these freedoms as long as you follow the license terms.

Under the following terms:

* Attribution - You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any 		  way that suggests the licensor endorses you or your use.
* NonCommercial - You may not use the material for commercial purposes.
* NoDerivatives - If you remix, transform, or build upon the material, you may not distribute the modified material.
* No additional restrictions - You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.
To view a full copy of this license, visit [this page](https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode).


Installation
============

## Requirements:

* Pyhon2.7: I recommend installing Anaconda Python Distribution from https://docs.continuum.io/anaconda/install that is platform independent.

Download the code (clone the repositoty) from github:

    git clone https://github.com/tevang/4D-CHAINS.git

Or to update an existing repository:
	
    git pull origin master

The next step is to install all Python dependencies. Enter the main directory and type:
	
    python setup.py install

or if you need sudo priviledges to install new python packages:
	
    sudo python setup.py install

Place the bin/ folder in your PATH environment variable if you want to call 4D-CHAINS executables from every directory. E.g. place the following line in your .bashrc:
	
    export PATH=<path to 4D-CHAINS/>/bin:$PATH

To run 4D-CHAINS use the 4Dchains.py script. You can do all operations you wish with this script as long as you provide the appropriate protocol file. E.g.
	
    4Dchains.py -protocol protocol.txt

You can generate a sample protocol file like this:
	
    4Dchains.py -writeprotocol

For information about all the available directives in the protocol file, please refer to the Manual.


Tutorials
============

You can find tutorials for Tudor (60 residues) and nEIt (248 residues) proteins under tutorials/ directory. In each folder you will find a protocol file and all the required input files. To run full protein assignment for the Tudor, do:

	4Dchains.py -protocol Tudor_protocol.txt

Likewise, for nEIt:

	4Dchains.py -protocol nEIt_protocol.txt

At the end you will find several output files and a directory named 4DCHAINS_workdir, which containts all the intermediate files for backtracking. To see what
these output files mean, do:

4Dchains.py -h






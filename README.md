

Brief Description
=================

4D-CHAINS is a software for fully automated protein backbone N-H and sidechain aliphatic C & H chemical shift assignment from 2 NMR spectra: a 4D TOCSY and
a 4D NOESY. Read our paper in [Nature Communications](https://www.nature.com/articles/s41467-017-02592-z) for more details.

<table style="border-collapse: collapse">
<tr>
<td style="vertical-align: top" valign="top">
    <strong>Abstract from the paper</strong>
    <p>Automated methods for NMR structure determination of proteins are continuously becoming more robust. However, current methods 
    addressing larger, more complex targets rely on analyzing 6–10 complementary spectra, suggesting the need for alternative 
    approaches. Here, we describe 4D-CHAINS/autoNOE-Rosetta, a complete pipeline for NOE-driven structure determination of medium- 
    to larger-sized proteins. The 4D-CHAINS algorithm analyzes two 4D spectra recorded using a single, fully protonated protein 
    sample in an iterative ansatz where common NOEs between different spin systems supplement conventional through-bond 
    connectivities to establish assignments of sidechain and backbone resonances at high levels of completeness and with a minimum 
    error rate. The 4D-CHAINS assignments are then used to guide automated assignment of long-range NOEs and structure refinement 
    in autoNOE-Rosetta. Our results on four targets ranging in size from 15.5 to 27.3 kDa illustrate that the structures of proteins 
    can be determined accurately and in an unsupervised manner in a matter of days.
    <p>
    <strong>Link to the paper</strong><br />
    <a href="https://www.nature.com/articles/s41467-017-02592-z">Nat. Commun</a>
    </p>
</td><td width="300">
<img src="images/network.png" width="300" /></img>
</td>
</tr>
</table>


License
============

4D-CHAINS software for protein NMR assignment is a property of **Thomas Evangelidis** and **Konstantinos Tripsianes** and is free for NON-COMMERCIAL USAGE. The code is licensed under the Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY- NC-ND 4.0). You are free to:

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

Place the bin/ folder in your PATH and lib/ in your PYTHONPATH environment variable if you want to call 4D-CHAINS executables from every directory. E.g. place the following line in your .bashrc:
	
    export PATH=<path to 4D-CHAINS/>/bin:$PATH
    export PYTHONPATH=<path to 4D-CHAINS/>/lib:$PYTHONPATH

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






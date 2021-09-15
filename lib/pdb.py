#!/usr/bin/env python
# 4D-CHAINS software is a property of is a property of Masaryk university and the authors are
# Thomas Evangelidis and Konstantinos Tripsianes. The code is licensed under the Attribution-NonCommercial-NoDerivatives 4.0
# International (CC BY-# NC-ND 4.0). You are free to:
# * Share - copy and redistribute the material in any medium or format.
# * The licensor cannot revoke these freedoms as long as you follow the license terms.
# Under the following terms:
# * Attribution - You must give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way 
#   that suggests the licensor endorses you or your use.
# * NonCommercial - You may not use the material for commercial purposes.
# * NoDerivatives - If you remix, transform, or build upon the material, you may not distribute the modified material.
# * No additional restrictions - You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.
# To view a full copy of this license, visit https://creativecommons.org/licenses/by-nc-nd/4.0/legalcode.


# import pymol
# pymol.finish_launching(['pymol', '-cq'])
from pymol import cmd
from pymol import stored
# execfile('/home/thomas/.pymolrc')
from .global_func import *
import sys


class PDB():

    helical_HN_distances = []
    sheet_HN_distances = []
    loop_HN_distances = []

    def __init__(self, PDB_FOLDER=".", PDB_PATTERN=".*\.pdb"):
        # Load all PDB files that match the pattern
        self.pdb_list = list_files(PDB_FOLDER, PDB_PATTERN, full_path=True)

    def get_chains(self, model):
        # find all the chain of this structure file
        stored.chains = []
        cmd.iterate(model, 'stored.chains.append(chain)')
        chains = list(set(stored.chains))
        chains.sort()
        return chains

    def save_helical_HN_distances(self, model, chain_sel):
        stored.resids = []
        cmd.select("helices", model + chain_sel + " and ss H and not resn PRO and name CA")
        cmd.iterate("helices", 'stored.resids.append(resi)')
        stored.resids = [int(r) for r in stored.resids]
        stored.resids.sort()
        # print "DEBUG: stored.resids=", stored.resids
        for i, resid in enumerate(stored.resids[0:-1]):
            next_resid = stored.resids[i + 1]
            if not next_resid == resid + 1:
                continue
            distNH = cmd.distance("distNH", model + chain_sel + "and resi "+str(resid)+" and name HN",
                                  model + chain_sel + " and resi "+str(next_resid)+" and name HN")
            self.helical_HN_distances.append(distNH)

    def save_sheet_HN_distances(self, model, chain_sel):
        stored.resids = []
        cmd.select("sheets", model + chain_sel + " and ss S and not resn PRO and name CA")
        cmd.iterate("sheets", 'stored.resids.append(resi)')
        stored.resids = [int(r) for r in stored.resids]
        stored.resids.sort()
        # print "DEBUG: stored.resids=", stored.resids
        for i, resid in enumerate(stored.resids[0:-1]):
            next_resid = stored.resids[i + 1]
            if not next_resid == resid + 1:
                continue
            distNH = cmd.distance("distNH", model + chain_sel + "and resi "+str(resid)+" and name HN",
                                  model + chain_sel + " and resi "+str(next_resid)+" and name HN")
            self.sheet_HN_distances.append(distNH)

    def save_loop_HN_distances(self, model, chain_sel):
        stored.resids = []
        cmd.select("lops", model + chain_sel + " and ss L and not resn PRO and name CA")
        cmd.iterate("lops", 'stored.resids.append(resi)')
        stored.resids = [int(r) for r in stored.resids]
        stored.resids.sort()
        # print "DEBUG: stored.resids=", stored.resids
        for i, resid in enumerate(stored.resids[0:-1]):
            next_resid = stored.resids[i + 1]
            if not next_resid == resid + 1:
                continue
            distNH = cmd.distance("distNH", model + chain_sel + "and resi "+str(resid)+" and name HN",
                                  model + chain_sel + " and resi "+str(next_resid)+" and name HN")
            self.loop_HN_distances.append(distNH)
    
    def measure_HN_distances(self):
        for pdb in self.pdb_list:
            model = os.path.basename(pdb).replace(".pdb", "")
            print("Loading file", pdb, "as", model)
            cmd.load(pdb, object=model)
            chains = self.get_chains(model)
            # print "DEBUG: chains=", chains
            for chain in chains:
                if chain:
                    chain_sel = " chain " + chain + " "
                else:
                    chain_sel = " "
                self.save_helical_HN_distances(model, chain_sel)
                self.save_sheet_HN_distances(model, chain_sel)
                self.save_loop_HN_distances(model, chain_sel)
            cmd.delete(model)

    def plot_HN_distance_distributions(self):
        self.helical_HN_distances.sort()
        self.sheet_HN_distances.sort()
        self.loop_HN_distances.sort()
        print("helical_HN_distances=", self.helical_HN_distances)
        print("sheet_HN_distances=", self.sheet_HN_distances)
        print("loop_HN_distances=", self.loop_HN_distances)
        helical_HN_distances = [d for d in helical_HN_distances if d!=0.0]  # remove 0 distances
        helical_hist, bin_edges = np.histogram(helical_HN_distances, bins=60, range=(0, 6.))
        helical_hist = helical_hist / float(helical_hist.sum())
        sheet_HN_distances = [d for d in sheet_HN_distances if d != 0.0]  # remove 0 distances
        sheet_hist, bin_edges = np.histogram(sheet_HN_distances, bins=60, range=(0, 6.))
        sheet_hist = sheet_hist / float(sheet_hist.sum())
        loop_HN_distances = [d for d in loop_HN_distances if d != 0.0]  # remove 0 distances
        loop_hist, bin_edges = np.histogram(loop_HN_distances, bins=60, range=(0, 6.))
        loop_hist = loop_hist / float(loop_hist.sum())

        import numpy
        import matplotlib.pyplot as plt
        import plotly.plotly as py  # tools to communicate with Plotly's server

        bins = numpy.linspace(0, 6, 60)   # 0.1 Angstrom bin width

        labels = ("SHEET", "HELIX", "LOOP")
        plt.bar(bins, sheet_hist, width=0.1, alpha=0.5)
        plt.bar(bins, helical_hist, width=0.1, alpha=0.5)
        plt.bar(bins, loop_hist, width=0.1, alpha=0.5)
        plt.xlabel('HN(i)-HN(i+1) distance in Angstroms', fontsize=24);
        plt.ylabel('Probability', fontsize=24);
        plt.legend(labels=labels)
        plt.show()

"""
# EXAMPLE:
# to run the code:
# pymol -cq pdb.py
pdb = PDB("/home/thomas/Downloads/talos/pdb")
pdb.measure_HN_distances()
pdb.plot_HN_distance_distributions()
sys.exit(1)
"""

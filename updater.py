import requests
import subprocess
import urllib.request
import time
import os
import logging
import sys
import argparse
import configparser

class UpdateDatabase(object):
    def main(self):
        """Main Program, updates the database"""
        print("Main")

    def getrmlsthelper(self, referencefilepath, update, start):
        """
        Makes a system call to rest_auth.pl, a Perl script modified from
        https://github.com/kjolley/BIGSdb/tree/develop/scripts/test
        And downloads the most up-to-date rMLST profile and alleles
        """
        from subprocess import call
        analysistype = 'rMLST'
        # Folders are named based on the download date e.g 2016-04-26
        # Find all folders (with the trailing / in the glob search) and remove the trailing /
        lastfolder = sorted(glob('{}{}/2*/'.format(referencefilepath, analysistype)))[-1].rstrip('/')
        delta, foldersize, d1 = schemedate(lastfolder)
        # Extract the path of the current script from the full path + file name
        homepath = os.path.split(os.path.abspath(__file__))[0]
        # Set the path/name of the folder to contain the new alleles and profile
        newfolder = '{}{}/{}'.format(referencefilepath, analysistype, d1)
        # System call
        rmlstupdatecall = 'cd {} && perl {}/rest_auth.pl -a {}/secret.txt'.format(newfolder, homepath, homepath)
        if update:
            if delta.days > 7 or foldersize < 100:
                printtime("Last update of rMLST profile and alleles was {} days ago. Updating".format(str(delta.days)),
                          start)
                # Create the path
                make_path(newfolder)
                # Copy over the access token to be used in the authentication
                shutil.copyfile('{}/access_token'.format(homepath), '{}/access_token'.format(newfolder))
                # Run rest_auth.pl
                call(rmlstupdatecall, shell=True)
                # Get the new alleles into a list, and create the combinedAlleles file
                alleles = glob('{}/*.tfa'.format(newfolder))
                combinealleles(start, newfolder, alleles)
            # If the profile and alleles are up-to-date, set :newfolder to :lastfolder
            else:
                newfolder = lastfolder
            # Ensure that the profile/alleles updated successfully
            # Calculate the size of the folder by adding the sizes of all the files within the folder together
            newfoldersize = sum(os.path.getsize('{}/{}'.format(newfolder, f)) for f in os.listdir(newfolder)
                                if os.path.isfile('{}/{}'.format(newfolder, f)))
            # If the profile/allele failed, remove the folder, and use the most recent update
            if newfoldersize < 100:
                shutil.rmtree(newfolder)
                newfolder = sorted(glob('{}{}/*/'.format(referencefilepath, analysistype)))[-1].rstrip('/')
        # Don't update the profile/alleles if not requested
        else:
            newfolder = lastfolder
        # Return the system call and the folder containing the profile and alleles
        return rmlstupdatecall, newfolder

    def getmlsthelper(self, referencefilepath, start, organism, update):
        """Prepares to run the getmlst.py script provided in SRST2"""
        from accessoryFunctions import GenObject
        # Initialise a set to for the organism(s) for which new alleles and profiles are desired
        organismset = set()
        # Allow for Shigella to use the Escherichia MLST profile/alleles
        organism = organism if organism != 'Shigella' else 'Escherichia'
        # As there are multiple profiles for certain organisms, this dictionary has the schemes I use as values
        organismdictionary = {'Escherichia': 'Escherichia coli#1',
                              'Shigella': 'Escherichia coli#1',
                              'Vibrio': 'Vibrio parahaemolyticus',
                              'Campylobacter': 'Campylobacter jejuni',
                              'Listeria': 'Listeria monocytogenes',
                              'Bacillus': 'Bacillus cereus',
                              'Klebsiella': 'Klebsiella pneumoniae'}
        # Allow for a genus not in the dictionary being specified
        try:
            organismset.add(organismdictionary[organism])
        except KeyError:
            # Add the organism to the set
            organismset.add(organism)
        for scheme in organismset:
            organismpath = os.path.join(referencefilepath, 'MLST', organism)
            # Find all folders (with the trailing / in the glob search) and remove the trailing /
            try:
                lastfolder = sorted(glob('{}/*/'.format(organismpath)))[-1].rstrip('/')
            except IndexError:
                lastfolder = []
            # Run the method to determine the most recent folder, and how recently it was updated
            delta, foldersize, d1 = schemedate(lastfolder)
            # Set the path/name of the folder to contain the new alleles and profile
            newfolder = '{}/{}'.format(organismpath, d1)
            if update:
                if delta.days > 7 or foldersize < 100:
                    printtime('Downloading {} MLST scheme from pubmlst.org'.format(organism), start)
                    # Create the object to store the argument attributes to feed to getmlst
                    getmlstargs = GenObject()
                    getmlstargs.species = scheme
                    getmlstargs.repository_url = 'http://pubmlst.org/data/dbases.xml'
                    getmlstargs.force_scheme_name = False
                    getmlstargs.path = newfolder
                    # Create the path to store the downloaded
                    make_path(getmlstargs.path)
                    getmlst.main(getmlstargs)
                    # Even if there is an issue contacting the database, files are created, however, they are populated
                    # with XML strings indicating that the download failed
                    # Read the first character in the file
                    try:
                        profilestart = open(glob('{}/*.txt'.format(newfolder))[0]).readline()
                    except IndexError:
                        profilestart = []
                    # If it is a <, then the download failed
                    if not profilestart or profilestart[0] == '<':
                        # Delete the folder, and use the previous definitions instead
                        shutil.rmtree(newfolder)
                        newfolder = lastfolder
                # If the profile and alleles are up-to-date, set :newfolder to :lastfolder
                else:
                    newfolder = lastfolder
            # If update isn't specified, don't update
            else:
                newfolder = lastfolder
                # Ensure that the profile/alleles updated successfully
                # Calculate the size of the folder by adding the sizes of all the files within the folder together
            try:
                newfoldersize = sum(os.path.getsize('{}/{}'.format(newfolder, f)) for f in os.listdir(newfolder)
                                    if os.path.isfile('{}/{}'.format(newfolder, f)))
            except (OSError, TypeError):
                newfoldersize = 100
            # If the profile/allele failed, remove the folder, and use the most recent update
            if newfoldersize < 100:
                shutil.rmtree(newfolder)
                try:
                    newfolder = sorted(glob('{}/*/'.format(organismpath)))[-1].rstrip('/')
                except IndexError:
                    newfolder = organismpath
            # Return the name/path of the allele-containing folder
            return newfolder

    def __init__(self):
        print("initialising")

if __name__ == '__main__':
    UpdateDatabase()

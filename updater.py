import requests
import subprocess
import urllib.request
import time
import os
import shutil
from Bio import SeqIO
import logging
import sys

import configparser
from glob import glob
import getmlst
from accessoryfunctions.accessoryFunctions import *
from pyaccessories.SaveLoad import SaveLoad

class UpdateDatabase(object):
    def main(self):
        """Main Program, updates the database"""
        print("Main")
        start_time = time.time()
        self.getrmlsthelper(start_time)

        for organism in self.loader.to_update:
            self.getmlsthelper(start_time, organism)

    def getrmlsthelper(self, start):
        """
        Makes a system call to rest_auth.pl, a Perl script modified from
        https://github.com/kjolley/BIGSdb/tree/develop/scripts/test
        And downloads the most up-to-date rMLST profile and alleles
        """
        from subprocess import call
        # Folders are named based on the download date e.g 2016-04-26
        # Find all folders (with the trailing / in the glob search) and remove the trailing /
        try:
            lastfolder = sorted(glob('{}{}/2*/'.format(self.referencefilepath, self.analysistype)))[-1].rstrip('/')
        except IndexError:
            lastfolder = "2000-01-01"

        delta, foldersize, d1 = self.schemedate(lastfolder)
        # Extract the path of the current script from the full path + file name
        homepath = os.path.split(os.path.abspath(__file__))[0]
        # Set the path/name of the folder to contain the new alleles and profile
        newfolder = '{}{}/{}'.format(self.referencefilepath, self.analysistype, d1)
        # System call
        rmlstupdatecall = 'cd {} && perl {}/rest_auth.pl -a {}/secret.txt'.format(newfolder, homepath, homepath)
        if foldersize < 100:
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
            self.combinealleles(start, newfolder, alleles)
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
            try:
                newfolder = sorted(glob('{}{}/*/'.format(self.referencefilepath, self.analysistype)))[-1].rstrip('/')
            except IndexError:
                pass
        # Return the system call and the folder containing the profile and alleles
        return rmlstupdatecall, newfolder

    def getmlsthelper(self, start, organism):
        """Prepares to run the getmlst.py script provided in SRST2"""
        # Initialise a set to for the organism(s) for which new alleles and profiles are desired
        organismset = set()
        # Allow for Shigella to use the Escherichia MLST profile/alleles
        organism = organism if organism != 'Shigella' else 'Escherichia'
        # As there are multiple profiles for certain organisms, this dictionary has the schemes I use as values

        # Allow for a genus not in the dictionary being specified
        try:
            organismset.add(self.loader.organismdictionary[organism])
        except KeyError:
            # Add the organism to the set
            organismset.add(organism)
        for scheme in organismset:
            organismpath = os.path.join(self.referencefilepath, 'MLST', organism)
            # Find all folders (with the trailing / in the glob search) and remove the trailing /
            try:
                lastfolder = sorted(glob('{}/*/'.format(organismpath)))[-1].rstrip('/')
            except IndexError:
                lastfolder = []
            # Run the method to determine the most recent folder, and how recently it was updated
            delta, foldersize, d1 = self.schemedate(lastfolder)
            # Set the path/name of the folder to contain the new alleles and profile
            newfolder = '{}/{}'.format(organismpath, d1)

            if foldersize < 100:
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

    @staticmethod
    def combinealleles(start, allelepath, alleles):
        printtime('Creating combined rMLST allele file', start)
        records = []

        # Open each allele file
        for allele in sorted(alleles):
            # with open(allele, 'rU') as fasta:
            for record in SeqIO.parse(open(allele, "rU"), "fasta"):
                # Extract the sequence record from each entry in the multifasta
                # Replace and dashes in the record.id with underscores
                record.id = record.id.replace('-', '_')
                # Remove and dashes or 'N's from the sequence data - makeblastdb can't handle sequences
                # with gaps
                # noinspection PyProtectedMember
                record.seq._data = record.seq._data.replace('-', '').replace('N', '')
                # Clear the name and description attributes of the record
                record.name = ''
                record.description = ''
                # Write each record to the combined file
                # SeqIO.write(record, combinedfile, 'fasta')
                records.append(record)
        with open('{}/rMLST_combined.fasta'.format(allelepath), 'w') as combinedfile:
            SeqIO.write(records, combinedfile, 'fasta')

    def schemedate(self, lastfolder):
        from datetime import date
        try:
            # Extract the folder name (date) from the path/name
            lastupdate = os.path.split(lastfolder)[1]
        except AttributeError:
            lastupdate = '2000-01-01'
        try:
            # Calculate the size of the folder by adding the sizes of all the files within the folder together
            foldersize = sum(os.path.getsize('{}/{}'.format(lastfolder, f)) for f in os.listdir(lastfolder)
                             if os.path.isfile('{}/{}'.format(lastfolder, f)))
        except (TypeError, FileNotFoundError):
            foldersize = 0
        # Try to figure out the year, month, and day from the folder name
        try:
            (year, month, day) = lastupdate.split("-")
            # Create a date object variable with the year, month, and day
            d0 = date(int(year), int(month), int(day))
        except ValueError:
            # Set an arbitrary date in the past to force an update
            d0 = date(2000, 1, 1)
        # Create a date object with the current date
        d1 = date(int(time.strftime("%Y")), int(time.strftime("%m")), int(time.strftime("%d")))
        # Subtract the last update date from the current date
        delta = d1 - d0

        return delta, foldersize, d1

    def __init__(self, parser):
        print("initialising")
        self.analysistype = "rMLST"
        # self.referencefilepath = "/mnt/nas/Adam/assemblypipeline/rMLST/"
        self.referencefilepath = os.path.join(parser.referencedirectory, "")
        self.start = parser.start
        self.loader = SaveLoad()

        # If the file was empty and it couldn't load but created the file
        import json
        try:
            # Fresh file
            if not self.loader.load("bacteria.json", True):
                self.loader.organismdictionary = {'Escherichia': 'Escherichia coli#1',
                                                  'Shigella': 'Escherichia coli#1',
                                                  'Vibrio': 'Vibrio parahaemolyticus',
                                                  'Campylobacter': 'Campylobacter jejuni',
                                                  'Listeria': 'Listeria monocytogenes',
                                                  'Bacillus': 'Bacillus cereus',
                                                  'Klebsiella': 'Klebsiella pneumoniae'}
                self.loader.to_update = list(self.loader.organismdictionary.keys())
                self.loader.dump("bacteria.json")

            if "organismdictionary" not in self.loader.__dict__:
                raise NameError

        except (json.decoder.JSONDecodeError, NameError):
            print("Invalid config file, please delete or fix")
            sys.exit(1)

        self.main()

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(description="descr")

    parser.add_argument('-r', '--referencedirectory',
                        required=True,
                        help='ref dir')
    parser.add_argument('-u', '--to_update',
                        required=False,
                        help='Comma separated list of bacteria to update')
    args = parser.parse_args()
    args.start = time.time()
    UpdateDatabase(args)

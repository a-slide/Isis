"""
@package    Utilities
@brief      Contains several usefull functions to interact with OS environement and to parse files
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~COMMAND LINE UTILITIES~~~~~~~#

def run_command(cmd, stdinput=None, control_stderr=True, ret_stderr=False, ret_stdout=True):
    """
    Run a command line in the default shell and return the standard output
    @param  cmd A command line string formated as a string
    @param  stdinput    Facultative parameters to redirect an object to the standard input
    @param  control_stderr   If True the standard error output will be verified
    @param  ret_stderr  If True the standard error output will be returned
    @param  ret_stdout  If True the standard output will be returned
    @note If ret_stderr and ret_stdout are True a tuple will be returned and if both are False
    None will be returned
    @return If no standard error return the standard output as a string
    @exception  OSError Raise if a message is return on the standard error output
    @exception  (ValueError,OSError) May be raise by Popen
    """
    # Function specific imports
    from subprocess import Popen, PIPE

    # Common error message in case of Exception
    msg = "An error occured while trying to execute the following command :\n{}".format(cmd)

    # Try to execute the command line in the default shell
    try:
        if stdinput:
            stdout, stderr = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE).communicate(input=stdinput)
        else:
            stdout, stderr = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()

        # If a message is returned on the std err output an SystemError is raised
        if control_stderr and stderr:
            raise OSError("{}\n{}".format(msg, stderr))

        # Else return data according to user choices is returned
        if ret_stdout and ret_stderr:
            return stdout, stderr
        elif ret_stdout:
            return stdout
        elif ret_stderr:
            return sterr
        else:
            return None

    # Take care of possible exceptions returned by Popen
    except (OSError, ValueError) as E:
        print (E)

def make_cmd_str(prog_name, opt_dict={}, opt_list=[]):
    """
    Create a Unix like command line string from a
    @param prog_name Name (if added to the system path) or path of the programm
    @param opt_dict Dictionnary of option arguments such as "-t 5". The option flag have to
    be the key (without "-") and the the option value in the dictionnary value. If no value is
    requested after the option flag "None" had to be asigned to the value field.
    @param opt_list List of simple command line arguments
    @exemple make_cmd_str("bwa", {"b":None, t":6, "i":"../idx/seq.fa"}, ["../read1", "../read2"])
    """

    # Start the string by the name of the program
    cmd = "{} ".format(prog_name)

    # Add options arguments from opt_dict
    if opt_dict:
        for key, value in opt_dict.items():
            if value:
                cmd += "-{} {} ".format(key, value)
            else:
                cmd += "-{} ".format(key)

    # Add arguments from opt_list
    if opt_list:
        for value in opt_list:
            cmd += "{} ".format(value)

    return cmd

#~~~~~~~FILE MANIPULATION~~~~~~~#

def mkdir(fp):
    """
    Create a directory at the indicated path\n
    Reproduce the ability of UNIX "mkdir -p" command
    (ie if the path already exits no exception will be raised).
    @param  fp path name where the folder should be created
    @exception  OSError Can be raise by os.mkdir
    """
    # Function specific imports
    from os import mkdir, path

    if path.exists(fp) and path.isdir(fp):
        print ("'{}' already exist in the current directory".format(fp))
    else:
        print ("Creating '{}' in the current directory".format(fp))
        mkdir(fp)

def file_basename (path):
    """
    Return the basename of a file without folder location and extension
    @param path Filepath as a string
    """
    return path.rpartition('/')[2].partition('.')[0]

def file_extension (path):
    """
    Return the extension of a file.
    @param path Filepath as a string
    """
    return path.rpartition(".")[2]

def file_name (path):
    """
    Return the complete name of a file with the extension but without folder location
    @param path Filepath as a string
    """
    return path.rpartition("/")[2]

def dir_name (path):
    """
    Return the complete path where is located the file without the file name
    @param path Filepath as a string
    """
    return path.rpartition("/")[0].rpartition("/")[2]

#~~~~~~~FASTA UTILITIES~~~~~~~#

def import_seq(filename, col_type="dict", seq_type="fasta"):
    """
    Import sequences from a fasta files in a list of biopython SeqRecord
    @param filename Valid path to a fasta file. Can contains several sequences and can be gzipped
    @param col_type Type of the collection where SeqReccord entries will be added ("list" or "dict").
    @param seq_type Type of the sequence file to parse (see Biopython seqIO for supported format)
    @return A list or a dictionnary containing all seqReccord objects from the fastq file
    @exception IOError  Raise if the path in invalid or unreadeable
    """
    # Require the Third party package Biopython
    from Bio import SeqIO
    import gzip

    # Try to open the file fist gz compressed and uncompressed
    try:

        # Verify if the type of the input sequence is valid
        seq_type = seq_type.lower()
        allowed_seq = ["fasta", "genbank", "gb", "fastq-illumina", "fastq-solexa" , "fastq",
        "fastq-sanger", "embl ", "abi ", "seqxml", "sff", "uniprot-xml"]
        assert seq_type in allowed_seq , "The input file format have to be in the following list : "+ ", ".join(allowed_seq)

        # Verify if the type of the output collection is valid
        col_type = col_type.lower()
        allowed_types = ["dict", "list"]
        assert col_type  in allowed_types, "The output collection type have to be in the following list : "+ ", ".join(allowed_types)

        # Open gzipped or uncompressed file
        if file_extension(filename) == "gz":
            print("\tUncompressing and extracting data")
            handle = gzip.open(filename, "r")
        else:
            print("\tExtracting data")
            handle = open(filename, "r")

        # Create the collection
        if col_type == "list":
            seq_col = list(SeqIO.parse(handle, seq_type))
        else:
            seq_col = SeqIO.to_dict(SeqIO.parse(handle, seq_type))

        # Close file, verify if the collection is filled and returned it
        handle.close()
        assert seq_col, 'The collection contains no SeqRecord after file parsing. Exit'
        return seq_col

    except IOError as E:
        print('CRITICAL ERROR. The fasta file ' + filename + ' is not readable. Exit')
        exit

    except AssertionError as E:
        print (E)
        exit

#~~~~~~~GRAPHICAL UTILIES~~~~~~~#

def fill_between_graph (X, Y, basename="out", img_type="png", title=None, xlabel=None, ylabel=None, baseline=0):
    """
    Trace a generic fill between graph with matplotlib pyplot
    @param X List of values for x axis
    @param Y List of values for y axis
    @param title Title of graph (facultative)
    @param xlabel Label for x axis (facultative)
    @param ylabel Label for y axis (facultative)
    @param basename Output basename of the image file (Default "out")
    @param img_type Type of the image file (Default "png")
    """

    # Require the Third party package Biopython
    from matplotlib import pyplot as plt

    # Create a figure object and adding details
    fig = plt.figure(figsize=(15, 10), dpi=100)
    if title:
        plt.title(title)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)

    # Plot an area representing the coverage depth
    plt.fill_between(X, Y, baseline, facecolor='green', alpha=0.5)

    # Tweak spacing to prevent clipping of ylabel
    plt.subplots_adjust(left=0.15)

    # Export figure to file
    try:
        fig.savefig(basename+"."+img_type, format = img_type)
    except ValueError as E:
        print (E)
        print ("Saving file as png")
        fig.savefig(basename+".png", format = "png")

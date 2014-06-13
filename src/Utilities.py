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

def run_shell(cmd):
    """
    Run a command line in the default shell and return the standard output
    @param  cmd a command line string formated as a string
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
        stdout, stderr = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()
        # If a message is returned on the std err output an SystemError is raised
        if stderr:
            raise OSError("{}\n{}".format(msg, stderr))
        # Else the std output data is returned
        else:
            return stdout

    # Take care of possible exceptions returned by Popen
    except (OSError, ValueError) as E:
        print (E)


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

    try:
        if path.exists(fp) and path.isdir(fp):
            print ("'{}' already exist in the current directory".format(fp))
        else:
            print ("Creating '{}' in the current directory".format(fp))
            mkdir(fp)

    # In case
    except OSError as E:
        print E
        raise

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

def import_fasta(filename, col_type="dict"):
    """
    Import sequences from a fasta files in a list of biopython SeqRecord
    @param filename Valid path to a fasta file. Can contains several sequences and can be gzipped
    @param col_type Type of the collection where SeqReccord entries will be added "list" or "dict".
    @return A list or a dictionnary containing all seqReccord objects from the fastq file
    @exception IOError  Raise if the path in invalid or unreadeable
    """

    # Require the Third party package Biopython
    from Bio import SeqIO
    import gzip

    # Try to open the file fist gz compressed and uncompressed
    try:
        if file_extension(filename) == "gz":
            print("\tUncompressing and extracting data")
            handle = gzip.open(filename, "r")
        else:
            print("\tExtracting data")
            handle = open(filename, "r")

        if col_type == "list":
            seq_col = list(SeqIO.parse(handle, "fasta"))
        else:
            seq_col = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

        handle.close()
        return seq_col

    # Try to open the file fist gz compressed and uncompressed
    except IOError as E:
        print('CRITICAL ERROR. The fasta file ' + filename + ' is not readable. Exit')
        exit
